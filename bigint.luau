--!strict

--[[
    Maintainers of Luau do not intend on supporting UInt64.
    https://github.com/luau-lang/luau/discussions/242#discussioncomment-12299931

    As a result, we are limited to using the provided bit32 library for
    bitwise operations, and using a base of 2^16 for our BigInt implementation.
]]--

local BigInt = {}
BigInt.__index = BigInt
type BigIntData = { Negative: boolean, Words: { number }}
export type BigInt = typeof(setmetatable({} :: BigIntData, BigInt))
export type Biggish = BigInt | number | string

-- BigInt algorithms
local function degenerate_divide (x: BigInt,
                                  y: number,
                                  should_return_remainder: boolean?
                                 ): BigInt
    -- Degenerate case of Algorithm D when "y" is not a BigInt
    if x:eq(0) then
        return BigInt:new(0)
    end

    local Q = table.create(#x.Words, 0)
    local k = 0

    for j = #Q, 1, -1 do
        local t = bit32.lshift(k, 16) + x.Words[j]
        Q[j] = t // y
        k = t - Q[j] * y
    end

    local result = BigInt:new()
    result.Negative = not (x.Negative == (y < 0))

    if should_return_remainder then
        result.Words = {
            bit32.band(k, 0xFFFF),
            bit32.band(bit32.rshift(k, 16), 0xFFFF)
        }
    else
        result.Words = Q
    end

    result:emsmallen()

    if result:eq(0) then
        result.Negative = false
    end

    return result
end

local function divide (x: BigInt,
                       y: BigInt,
                       should_return_remainder: boolean?
                      ): BigInt
    --[[
        Algorithm D from The Art of Computer Programming
        Vol. 2 Seminumerical Algorithms, 3rd ed. pg. 272

        See also: Hacker's Delight divmnu.c
    ]]--

    -- Assumption is that x and y are base 2^16 digits of nonnegative numbers.

    -- We don't return the remainder with every division call as it requires a
    -- nontrivial denormalization operation which incurs additonal computation.

    if x:eq(0) then
        return BigInt:new(0)
    end

    if (#y.Words == 1) then
        if (y.Words[1] == 0) then
            error("Div by zero")
        end

        return degenerate_divide (x, y.Words[1], should_return_remainder)
    end

    if x:abs():__lt(y:abs()) then
        if should_return_remainder then
            return x
        end

        return BigInt:new(0)
    end

    if x:abs():__eq(y:abs()) then
        if should_return_remainder then
            return BigInt:new(0)
        end

        return BigInt:new(1)
    end

    local Q = table.create(#x.Words, 0)
    local S = bit32.countlz(y.Words[#y.Words]) - 16
    local U = x:shl(S).Words
    local V = y:shl(S).Words

    table.insert(U, 0) -- Ensure #U > #V
    for j = (#U - #V), 1, -1 do
        local temp = bit32.bor(bit32.lshift(U[#V + j], 16), U[#V + j - 1])
        -- TODO: Check whether using `rhat := temp % V[#V]` instead of
        --       what we have now that the Luau JIT will be able to optimize
        --       this down into one div instruction.
        local qhat = temp // V[#V]
        local rhat = temp - qhat * V[#V]

        repeat
            local cond1 = qhat >= 2 ^ 16
            local upperbound = bit32.bor(bit32.lshift(rhat, 16), U[#V + j - 2])
            local cond2 = (qhat * V[#V - 1]) > upperbound

            if cond1 or cond2 then
                qhat -= 1
                rhat += V[#V]

                continue
            end

            break
        until rhat >= 2 ^ 16

        local k = 0
        local t = 0

        for i = 1, #V do
            local p = qhat * V[i]
            t = U[i + j - 1] - bit32.band(p, 0xFFFF) - k
            U[i + j - 1] = bit32.band(t, 0xFFFF)
            k = bit32.band(bit32.rshift(p, 16) - bit32.rshift(t, 16), 0xFFFF)
        end

        t = U[#V + j] - k
        U[#V + j] = bit32.band(t, 0xFFFF)

        Q[j] = bit32.band(qhat, 0xFFFF)
        if t < 0 then
            Q[j] -= 1
            k = 0

            for i = 1, #V do
                t = U[i + j - 1] + V[i] + k
                U[i + j - 1] = bit32.band(t, 0xFFFF)
                k = bit32.rshift(t, 16)
            end

            U[#V + j] = bit32.band(U[#V + j] + k, 0xFFFF)
        end
    end

    local result = BigInt:new()
    result.Negative = not (x.Negative == y.Negative)

    if should_return_remainder then
        result.Words = U
        result:emsmallen()

        return result:shr(S)
    end

    result.Words = Q
    result:emsmallen()

    if result:eq(0) then
        result.Negative = false
    end

    return result
end

local function __divmod(div: (number) -> number,
                        mod: (number) -> number
                       ): (number) -> (number, number)
    return function (product: number): (number, number)
        return div(product), mod(product)
    end
end

local __divmod16 = __divmod(
    function (x) return bit32.rshift(x, 16) end,
    function (x) return bit32.band(x, 0xFFFF) end
)

local __divmod10 = __divmod(
    function (x) return x // 10000 end,
    function (x) return x % 10000 end
)

local function naive_muladd (x: {number},
                             mul: {number},
                             add: {number}?,
                             operation: ((number) -> (number, number))?
                            ): {number}
    --[[
        Algorithm M from The Art of Computer Programming
        Vol. 2 Seminumerical Algorithms, 3rd ed. pg. 268
    ]]--

    local result = {}
    local divmod = operation or __divmod16

    if add then
        result = table.clone(add)
    end

    -- Make result have size = (#x.Words + #mul.Words)
    local padding = #x + #mul - #result
    for _ = 1, padding do
        table.insert(result, 0)
    end

    for i = 1, #x do
        local carry = 0

        for j = 1, #mul do
            local product = x[i] * mul[j] + result[i + j - 1] + carry
            carry, result[i + j - 1] = divmod(product)
        end

        local temp = result[i + #mul] + carry
        result[i + #mul + 1], result[i + #mul] = divmod(temp)
    end

    for i = #result, 1, -1 do
        if (result[i] ~= 0) then break end

        if #result > 1 then
            table.remove(result, i)
        end
    end

    return result
end

function karatsuba (x: BigInt, y: BigInt): BigInt
    -- When using Karatsuba's method for numbers >= 2^1024 we see around a
    -- 30% performance improvement when multiplying two numbers of around equal
    -- size in bits.
    local m = math.max(#x, #y) // 2
    local x0 = x:take(m)
    local x1 = x:shr(m)
    local y0 = y:take(m)
    local y1 = y:shr(m)

    local z2 = x1 * y1
    local z0 = x0 * y0
    local z1 = (x1 + x0) * (y0 + y1) - z2 - z0

    return z2:shl(2 * m) + z1:shl(m) + z0
end

function BigInt.take (num: BigInt, bits: number): BigInt
    -- Take last N bits of a BigInt
    local result = num:clone()

    if #result <= bits then
        return result
    end

    local groups = (bits // 16) + 1
    local remainder = bits % 16

    if remainder > 0 then
        groups += 1
    end

    for i = #result.Words, groups, -1 do
        table.remove(result.Words, i)
    end

    if remainder > 0 then
        local temp = result.Words[#result.Words]
        local mask = (2 ^ remainder) - 1
        result.Words[#result.Words] = bit32.band(temp, mask)
    end

    return result
end

-- BigInt implementation
function BigInt:new (value: Biggish?,
                     radix: number?
                    ): BigInt
    if type(value) == "table" and getmetatable(value) == self then
        return value
    end

    local instance = setmetatable({} :: BigIntData, self)

    instance.Words = { 0 }
    instance.Negative = false

    if type(value) == "string" then
        return self.parse(value, radix)
    end

    if type(value) == "number" then
        instance.Negative = value < 0
        local magnitude = math.abs(value)

        if magnitude < 2 ^ 16 then
            instance.Words[1] = magnitude
            return instance
        end

        if magnitude >= 2 ^ 53 then
            error("BigInt:new() value passed as number literal exceeds 2 ^ 53")
        end

        table.remove(instance.Words)

        if magnitude < 2 ^ 32 then
            repeat
                table.insert(instance.Words, bit32.band(magnitude, 0xFFFF))
                magnitude = bit32.rshift(magnitude, 16)
            until magnitude == 0
        else
            repeat
                table.insert(instance.Words, magnitude % (2 ^ 16))
                magnitude = magnitude // (2 ^ 16)
            until magnitude == 0
        end
    end

    return instance
end

function BigInt:copy (other)
    self.Words = table.clone(other.Words)
    self.Negative = other.Negative
end

function BigInt.clone (num: BigInt): BigInt
    local result = BigInt:new()

    result:copy(num)

    return result
end

function BigInt.parse_hex (str: string): BigInt
    local negative = false
    str = str:gsub("%s", "")

    if str:match("^[+-]") ~= nil then
        negative = str:sub(1, 1) == "-"
        str = str:sub(2)
    end

    local num_words = #str // 4
    if #str % 4 > 0 then
        num_words += 1
    end

    local words = table.create(num_words, 0)
    local offset = 0

    if #str % 4 > 0 then
        local word = tonumber(str:sub(1, #str % 4), 16)

        if not word then
            error("Invalid value passed to BigInt.parse_hex()")
        end

        words[#words] = word
        offset = 1
        str = str:sub(#str % 4 + 1)
    end

    for i = 0, num_words - offset - 1 do
        local temp = str:sub(i * 4 + 1, (i + 1) * 4)
        local word = tonumber(temp, 16)

        if not word then
            error("Invalid value passed to BigInt.parse_hex()")
        end

        words[#words - i - offset] = word
    end

    local result = BigInt:new()
    result.Words = words
    result.Negative = negative

    return result
end

function BigInt.parse (str: string, radix: number?): BigInt
    -- TODO: Handle other radii

    if (str:match("^[+-]?0x") ~= nil) then
        return BigInt.parse_hex(str:sub(3))
    end

    if (radix == 16) then
        return BigInt.parse_hex(str)
    elseif radix and (radix ~= 10) then
        error("Unsupported radix " .. tostring(radix))
    end

    local negative = false
    if str:match("^[+-]") ~= nil then
        negative = str:sub(1, 1) == "-"
        str = str:sub(2)
    end

    str = str:gsub("[%s_',;.]", "")

    -- digits = floor(log_10(base)) = floor(log_10(2^16)) = 4
    local digits = 4
    local words = {}

    for i = #str, 1, -digits do
        local num = 0
        local place = 1

        local stop = math.max(1, i - digits + 1)
        for j = i, stop, -1 do
            local char = str:sub(j, j)

            -- TODO: handle and validate digits for other radii
            local digit = tonumber(char)

            if not digit then
                error("Invalid digit")
            end

            num += digit * place
            place *= 10
        end

        table.insert(words, num)
    end

    local z = BigInt:new()
    local base = BigInt:new(10 ^ digits)

    for i = #words, 1, -1 do
        local add = BigInt:new(words[i])
        z.Words = naive_muladd(z.Words, base.Words, add.Words)
    end

    z.Negative = negative

    return z
end

function BigInt.__tostring (num: BigInt)
    local result = ""

    if #num.Words == 1 then
        if num.Negative then
            result = "-" .. result
        end

        return result .. tostring(num.Words[1])
    end

    local z = BigInt:new()

    for i = #num.Words, 1, -1 do
        z.Words = naive_muladd(z.Words,
                               { 2 ^ 16 },
                               { num.Words[i] },
                               __divmod10)
    end

    for i = #z.Words, 1, -1 do
        local word = bit32.band(z.Words[i], 0xFFFF)
        result ..= string.format("%04d", word)
    end

    -- Replace leading zeros with nothing
    result = result:sub(1, #result - 1):gsub("^0+", "")
          .. result:sub(#result, #result)

    if num.Negative then
        result = "-" .. result
    end

    return result
end

function BigInt.hex (num: BigInt, length: number?): string
    local result = ""
    for i = #num.Words, 1, -1 do
        result ..= string.format("%04x", num.Words[i])
    end

    -- Replace leading zeros with nothing
    result = result:sub(1, #result - 1):gsub("^0+", "")
          .. result:sub(#result, #result)

    if length and length > #result then
        for i = 1, length - #result do
            result = "0" .. result
        end
    end

    return result
end

function BigInt:embiggen (other)
    if #self.Words > #other.Words then return end

    local padding = #other.Words - #self.Words + 1
    for i = 1, padding do
        table.insert(self.Words, 0)
    end
end

function BigInt:emsmallen ()
    if #self.Words == 1 then return end

    for i = #self.Words, 2, -1 do
        if (self.Words[i] ~= 0) then return end

        table.remove(self.Words, i)
    end

end

function BigInt.__len (num: BigInt): number
    return num:bits()
end

function BigInt.bits (num: BigInt): number
    -- Computes floor(log_2(num)) + 1, or the number of binary digits.
    -- TODO: We don't need to compute this every time we call this function.
    return ((#num.Words - 1) * 16) + (32 - bit32.countlz(num.Words[#num.Words]))
end

function BigInt.abs (num: BigInt): BigInt
    local result = num:clone()
    result.Negative = false

    return result
end

function BigInt.value (num: BigInt): number?
    if #num < 53 then
        local result = 0
        local exponent = 0

        for i = 1, #num.Words do
            result += num.Words[i] * (2 ^ exponent)
            exponent += 16
        end

        if num.Negative then
            return -result
        end

        return result
    end

    return nil
end

function BigInt.bit (num: BigInt, index: number): number?
    -- FIXME: Does it make more sense to return a boolean?
    if (index == 0) or (index > #num) then
        return nil
    end

    index -= 1
    local word = num.Words[index // 16 + 1]
    return bit32.band(bit32.rshift(word, index % 16), 1)
end

function BigInt.countrz (num: BigInt): number
    local zeros = 0
    for _, word in ipairs(num.Words) do
        if word == 0 then
            zeros += 16
        else
            zeros += bit32.countrz(bit32.band(word, 0xFFFF))
            break
        end
    end

    return zeros
end

function BigInt.is_power_of_2 (num: BigInt): boolean
    return (#num - num:countrz()) == 1
end

function BigInt.divide (_lhs: Biggish,
                        _rhs: Biggish,
                        should_return_remainder: boolean?
                       ): BigInt
    local lhs: BigInt = BigInt:new(_lhs)

    if (type(_rhs) == "number") and (math.abs(_rhs) < 2 ^ 16) then
        local result = degenerate_divide(lhs,
                                         math.abs(_rhs),
                                         should_return_remainder)
        local is_rhs_negative = math.sign(_rhs) == -1

        result.Negative = not (result.Negative == is_rhs_negative)
        return result
    end

    local rhs: BigInt = BigInt:new(_rhs)
    return divide(lhs, rhs, should_return_remainder)
end

function BigInt.__div (lhs: Biggish, rhs: Biggish): BigInt
    return BigInt.divide(lhs, rhs)
end

function BigInt.__idiv (lhs: Biggish, rhs: Biggish): BigInt
    return BigInt.divide(lhs, rhs)
end

function BigInt.__mod (lhs: Biggish, _rhs: Biggish): BigInt
    local rhs: BigInt = BigInt:new(_rhs)
    local result = BigInt.divide(lhs, rhs, true)

    if result.Negative ~= rhs.Negative then
        return result + rhs
    end

    return result
end

function BigInt.__mul (_lhs: Biggish, _rhs: Biggish): BigInt
    local lhs: BigInt = BigInt:new(_lhs)
    local rhs: BigInt = BigInt:new(_rhs)

    if lhs:eq(0) or rhs:eq(0) then
        return BigInt:new()
    end

    if lhs:eq(1) then
        return rhs:clone()
    end

    if rhs:eq(1) then
        return lhs:clone()
    end

    local result
    if rhs:is_power_of_2() then
        result = lhs:shl(#rhs - 1)
    elseif lhs:is_power_of_2() then
        result = rhs:shl(#lhs - 1)
    else
        -- TODO: Implement Toom-k methods of multiplication
        if ((#lhs >= 768) and (#rhs >= 512)) or
           ((#lhs >= 512) and (#rhs >= 768)) then
            result = karatsuba(lhs, rhs)
        else
            result = lhs:clone()
            result.Words = naive_muladd(result.Words, rhs.Words)
        end
    end

    result.Negative = not (lhs.Negative == rhs.Negative)
    result:emsmallen()

    return result
end

function BigInt.__pow (lhs: BigInt, _rhs: Biggish): BigInt
    -- Using fast powering algorithm.
    -- This uses the naive fast powering algorithm which is susceptible to
    -- timing side channel attacks.

    -- The exponent is limited to a range [0, 2 ^ 32] as it is improbable that
    -- one would attempt to raise something to a power above the limit.
    local rhs: BigInt = BigInt:new(_rhs)
    local exponent: number? = rhs:value()

    if (exponent == nil) or (exponent >= 2 ^ 32) or (exponent < 0) then
        error("Invalid exponent")
    end

    if exponent == 0 then
        return BigInt:new(1)
    end

    if exponent == 2 then
        return lhs * lhs
    end

    local x = lhs:clone()
    local y = BigInt:new(1)
    while exponent > 1 do
        if exponent % 2 == 1 then
            y = x * y
            exponent -= 1
        end

        x = x * x
        exponent = bit32.rshift(exponent, 1)
    end

    return x * y
end

function __add (u: {number}, v: {number}): {number}
    --- Algorithm A
    local k = 0

    local m = math.max(#u, #v) + 1
    local w = table.create(m, 0)

    for j = 1, m - 1 do
        local tmp = (u[j] or 0) + (v[j] or 0) + k

        w[j] = bit32.band(tmp, 0xFFFF)
        k = bit32.rshift(tmp, 16)
    end

    w[m] = k

    return w
end

function BigInt.__add (_lhs: Biggish, _rhs: Biggish): BigInt
    local lhs: BigInt = BigInt:new(_lhs)
    local rhs: BigInt = BigInt:new(_rhs)

    if rhs.Negative then
        return lhs - (-rhs)
    end

    if lhs.Negative then
        return -((-lhs) - rhs)
    end

    local result = lhs:clone()
    result.Words = __add(lhs.Words, rhs.Words)

    result:emsmallen()

    if result:eq(0) then
        result.Negative = false
    end

    return result
end

function BigInt.__sub (_lhs: Biggish, _rhs: Biggish): BigInt
    local lhs: BigInt = BigInt:new(_lhs)
    local rhs: BigInt = BigInt:new(_rhs)

    if rhs.Negative then
        return lhs + -rhs
    end

    if lhs.Negative then
        return -(lhs:abs() + rhs)
    end

    local result = BigInt:new()

    local lwords = lhs.Words
    local rwords = rhs.Words

    if lhs < rhs then
        lwords, rwords = rwords, lwords
        result.Negative = true

        result:embiggen(rhs)
    end

    local borrow = 0
    local i = 1
    local r = 1

    for l = 1, #lwords do
        local difference = lwords[l] + borrow

        if r ~= #rwords + 1 then
            difference -= rwords[r]
            r += 1
        end

        result.Words[i] = bit32.band(difference, 0xFFFF)
        borrow = bit32.arshift(difference, 16)

        i += 1
    end

    result:emsmallen()

    if result:eq(0) then
        result.Negative = false
    end

    return result
end

function BigInt.__unm (num: BigInt): BigInt
    local result = num:clone()
    result.Negative = not result.Negative

    return result
end

function BigInt.cmp (lhs: BigInt, rhs: BigInt): number
    -- Returns {-1, 0, 1} for gt, eq, lt like C++ <=> operator
    if lhs == rhs then
        return 0
    end

    if (lhs.Negative) and (not rhs.Negative) then
        return 1
    end

    if (not lhs.Negative) and (rhs.Negative) then
        return -1
    end

    if #lhs.Words > #rhs.Words then
        if lhs.Negative and rhs.Negative then
            return -1
        else
            return 1
        end
    end

    if #lhs.Words < #rhs.Words then
        if lhs.Negative and rhs.Negative then
            return 1
        else
            return -1
        end
    end

    if #lhs > #rhs then
        if lhs.Negative and rhs.Negative then
            return -1
        else
            return 1
        end
    end

    if #lhs < #rhs then
        if lhs.Negative and rhs.Negative then
            return 1
        else
            return -1
        end
    end

    for i = #lhs.Words, 1, -1 do
        local lword = lhs.Words[i]
        local rword = rhs.Words[i]

        if lword == rword then
            continue
        end

        if lword < rword then
            if lhs.Negative and rhs.Negative then
                return 1
            else
                return -1
            end
        end

        if lword > rword then
            if lhs.Negative and rhs.Negative then
                return -1
            else
                return 1
            end
        end
    end

    return 0
end

function BigInt.eq (_lhs: Biggish, _rhs: Biggish): boolean
    local lhs = BigInt:new(_lhs)

    if (type(_rhs) == "number") and (math.abs(_rhs) < 2 ^ 16) then
        if #lhs.Words > 1 then
            return false
        end

        if (math.sign(_rhs) == -1) and (not lhs.Negative) then
            return false
        end

        if (math.sign(_rhs) >= 0) and lhs.Negative then
            return false
        end

        return lhs.Words[1] == math.abs(_rhs)
    end

    local rhs: BigInt = BigInt:new(_rhs)
    return lhs:__eq(rhs)
end

function BigInt.neq (_lhs: Biggish, _rhs: Biggish): boolean
    return not BigInt.eq(_lhs, _rhs)
end

function BigInt.__eq (lhs: BigInt, rhs: BigInt): boolean
    if lhs.Negative ~= rhs.Negative then
        return false
    end

    if #lhs.Words ~= #rhs.Words then
        return false
    end

    if #lhs ~= #rhs then
        return false
    end

    for i = 1, #lhs.Words do
        if lhs.Words[i] ~= rhs.Words[i] then
            return false
        end
    end

    return true
end

function BigInt.lt (_lhs: Biggish, _rhs: Biggish): boolean
    local lhs = BigInt:new(_lhs)

    if (type(_rhs) == "number") and (math.abs(_rhs) < 2 ^ 16) then
        if #lhs.Words ~= 1 then
            return false
        end

        if lhs.Negative and math.sign(_rhs) >= 0 then
            return true
        end

        if not lhs.Negative and math.sign(_rhs) == -1 then
            return false
        end

        if lhs.Negative and math.sign(_rhs) == -1 then
            return lhs.Words[1] > _rhs
        end

        return lhs.Words[1] < _rhs
    end

    return lhs:__lt(BigInt:new(_rhs))
end

function BigInt.__lt (lhs: BigInt, rhs: BigInt): boolean
    return lhs:cmp(rhs) < 0
end

function BigInt.le (_lhs: Biggish, _rhs: Biggish): boolean
    return BigInt.eq(_lhs, _rhs) or BigInt.lt(_lhs, _rhs)
end

function BigInt.__le (lhs: BigInt, rhs: BigInt): boolean
    return lhs:cmp(rhs) <= 0
end

function BigInt.gt (_lhs: Biggish, _rhs: Biggish): boolean
    return not BigInt.le(_lhs, _rhs)
end

function BigInt.ge (_lhs: Biggish, _rhs: Biggish): boolean
    return not BigInt.lt(_lhs, _rhs)
end

-- Luau does not implement Lua 5.3 metamethods for bitwise operators
function BigInt.shl (lhs: BigInt, rhs: (BigInt | number)): BigInt
    local shift = 0
    if type(rhs) == "number" then
        shift = rhs
    else
        -- It's unlikely that there's a use case for shifting something more
        -- than 2^53 bits and we will most likely just encounter OOM errors.
        local value = rhs:value()

        if value then
            shift = value
        else
            error("Invalid shift")
        end
    end

    if shift == 0 then
        return lhs:clone()
    end

    if shift < 0 then
        return lhs:shr(-shift)
    end

    local result = lhs:clone()
    local groups = shift // 16
    for i = 1, groups do
        table.insert(result.Words, 1, 0)
    end

    local bit_shift = shift % 16
    local leading_zeros = bit32.countlz(lhs.Words[#lhs.Words]) - 16

    if leading_zeros < bit_shift then
        table.insert(result.Words, 0)
    end

    for i = #result.Words, 2, -1 do
        local shift_cur = bit32.lshift(result.Words[i], bit_shift)
        local shift_prev = bit32.rshift(result.Words[i - 1], 16 - bit_shift)
        local shifted = bit32.bor(shift_cur, shift_prev)

        result.Words[i] = bit32.band(shifted, 0xFFFF)
    end

    local shifted = bit32.lshift(result.Words[1], bit_shift)
    result.Words[1] = bit32.band(shifted, 0xFFFF)

    result:emsmallen()

    return result
end

function BigInt.shr (lhs: BigInt, rhs: (BigInt | number)): BigInt
    local shift = 0
    if type(rhs) == "number" then
        shift = rhs
    else
        -- It's highly unlikely there's a use case in which we shift anything
        -- more than 2 ^ 53 bits.
        local value = rhs:value()

        if value then
            shift = value
        else
            error("Invalid shift")
        end
    end

    if shift == 0 then
        return lhs:clone()
    end

    if shift < 0 then
        return lhs:shl(-shift)
    end

    if shift > #lhs then
        return BigInt:new(0)
    end

    local result = lhs:clone()
    local groups = shift // 16

    for i = 1, groups do
        table.remove(result.Words, 1)
    end

    local bit_shift = shift % 16
    for i = 1, #result.Words - 1 do
        local shift_cur = bit32.rshift(result.Words[i], bit_shift)
        local shift_next = bit32.lshift(result.Words[i + 1], 16 - bit_shift)
        local shifted = bit32.bor(shift_cur, shift_next)

        result.Words[i] = bit32.band(shifted, 0xFFFF)
    end

    local shifted = bit32.rshift(result.Words[#result.Words], bit_shift)
    result.Words[#result.Words] = bit32.band(shifted, 0xFFFF)
    result:emsmallen()

    return result
end

return BigInt
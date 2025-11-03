--!strict

local BigInt = require("./bigint")
type BigInt = BigInt.BigInt
type Biggish = BigInt.Biggish

local module = {
    internal = {},
    nt = { -- number theory
        internal = {}
    },
    prime = {
        internal = {}
    },
    random = {}
}

function module.random.random (bits: number?): BigInt
    -- Generates a random number from [0, 2 ^ bits - 1]
    -- If bits is nil, then generates a random 128-bit integer.
    -- The randomness depends on the underlying PRNG provided by math.random()
    -- in Luau, and should not be considered cryptographically secure.

    local result = BigInt:new()

    local temp = bits or 128
    local groups = temp // 16
    local remainder = temp % 16

    result.Words = table.create(groups)

    for i = 1, groups do
        result.Words[i] = math.random(2 ^ 16) - 1
    end

    if remainder > 0 then
        local value = math.random(2 ^ remainder) - 1
        table.insert(result.Words, value)
    end

    return result
end

function module.random.random_prime (bits: number?): BigInt
    -- Generate (n - 1) bits and append a 1 to generate an n-bit odd prime
    -- FIXME: Handle case for n=1 or n=2.
    local n: number = bits or 128

    while true do
        local num = module.random.random(n - 1):shl(1)
        if num:eq(0) then continue end

        num.Words[1] += 1

        local np = num - 1
        local r = np:countrz()
        local d = np:shr(r)
        if module.prime.internal.miller_rabin(num, r, d, 2, np) then
            -- MR with 2 as a witness is the fastest MR test we can perform.
            -- There is at most a 1 in 4 chance a composite number passes this
            -- step, so we can save our more computationally expensive checks
            -- for the numbers that pass this step.
            --
            -- This alone gave a ~35x speed up for generating 512-bit primes.
            continue
        end

        if module.prime.is_prime(num, 2 ^ -8) then
            return num
        end
    end
end


function module.prime.internal.p_kt (k: number, t: number, M: number)
    --[[
    According to FIPS 186-5 (doi:10.6028/NIST.FIPS.186-5)

    p_{k,t} is the probability that an odd k-bit integer which passes t rounds
    of Miller-Rabin is actually composite.

    Given an integer with k-bits and a p_target, we can compute a t such that
    p_{k,t} <= p_target.

    This function computes p_{k,t} for a given k, t, and M according to
    eq. (2) in pg. 72 of FIPS 186-5.
    ]]--
    local coefficient = 2.00743 * math.log(2) * k * (2 ^ -k)
    local sum_coefficient = 8 * ((math.pi ^ 2) - 6) / 3 * 2 ^ (k - 2)

    local sum = 0

    for m = 3, M do
        local w = 2 ^ (m - (m - 1) * t)

        local inner_sum = 0
        for j = 2, m do
            inner_sum += 2 ^ -(j + (k - 1) / j)
        end

        sum += w * inner_sum
    end

    sum = 2 ^ (k - 2 - M * t) + sum_coefficient * sum

    return coefficient * sum
end

function module.prime.internal.get_num_tests (k: number,
                                              p_target: number
                                             ): number
    -- We lose precision around k = 1000, so just hardcode 2 tests
    if k >= 1000 then return 2 end

    local T = math.ceil(-math.log(p_target, 2) / 2)
    local K = math.floor(2 * math.sqrt(k - 1) - 1)

    -- NOTE: This deviates from FIPS 186-5 in that t will never equal 1.
    --       This is because in the MR-primality test, we always perform an
    --       MR(2) test then choose subsequent tests randomly.
    --       While picking MR(2) is should be indistinguishable from randomly
    --       picking another witness, it gives me peace of mind to do another.
    --       Additionally, it is possible to adversarially generate probable
    --       primes such that they pass MR(2).

    for t = 2, T do
        for M = 3, K do
            if module.prime.internal.p_kt(k, t, M) <= p_target then
                return t
            end
        end
    end

    return T
end

-- List of primes to potentially be used as witnesses in Miller-Rabin
module.prime.internal.primes = {
      2,    3,    5,    7,   11,   13,   17,   19,   23,   29,   31,   37,   41,
     43,   47,   53,   59,   61,   67,   71,   73,   79,   83,   89,   97,  101,
    103,  107,  109,  113,  127,  131,  137,  139,  149,  151,  157,  163,  167,
    173,  179,  181,  191,  193,  197,  199,  211,  223,  227,  229,  233,  239,
    241,  251,  257,  263,  269,  271,  277,  281,  283,  293,  307,  311,  313,
    317,  331,  337,  347,  349,  353,  359,  367,  373,  379,  383,  389,  397,
    401,  409,  419,  421,  431,  433,  439,  443,  449,  457,  461,  463,  467,
    479,  487,  491,  499,  503,  509,  521,  523,  541,  547,  557,  563,  569,
    571,  577,  587,  593,  599,  601,  607,  613,  617,  619,  631,  641,  643,
    647,  653,  659,  661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
    739,  743,  751,  757,  761,  769,  773,  787,  797,  809,  811,  821,  823,
    827,  829,  839,  853,  857,  859,  863,  877,  881,  883,  887,  907,  911,
    919,  929,  937,  941,  947,  953,  967,  971,  977,  983,  991,  997, 1009,
   1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091,
   1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187,
   1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283,
   1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381,
   1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481,
   1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567,
   1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
   1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753,
   1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871,
   1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979,
   1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029 }

function module.prime.internal.miller_rabin (num: BigInt,
                                             r: number,
                                             d: BigInt,
                                             witness: Biggish,
                                             np: BigInt?
                                            ): boolean
    -- Test whether num is _composite_ against a given witness via Miller-Rabin
    -- Let r, d be defined such that num := d * 2 ^ r

    if num:eq(witness) then return true end
    local x = module.modexp(witness, d, num)

    for i = 1, r do
        local y = (x * x) % num

        if y:eq(1) and (not x:eq(1)) and (not x:eq(np or num - 1)) then
            return true
        end

        x = y
    end

    if not x:eq(1) then
        return true
    end

    return false -- not composite according to witness
end

function module.internal.sqrt_optimization (num: BigInt): BigInt
    local x = num:shr(#num // 2)

    -- Newton's method converges in about O(log N) so I arbitrarily chose
    -- a maximum of 2 log_2(N) steps for this test.
    for i = 1, (2 * #num) do
        x = (((x * x) + num) / x):shr(1)

        local distance = (num - (x ^ 2)):abs()
        if #distance.Words == 1 and distance.Words[1] <= 5 then
            -- https://math.stackexchange.com/questions/41337#comment91968_41355
            -- W.C. Jagy states that GMP uses this as the stop criterion.

            if distance.Words[1] == 0 then
                return x
            end

            break
        end
    end

    return x
end

function module.modsqrt (_num: Biggish, _p: Biggish): BigInt?
    --- Computes r such that its quadradic residue is n under mod p.
    local num: BigInt = BigInt:new(_num)
    local p: BigInt = BigInt:new(_p)

    if module.nt.jacobi(num, p) ~= 1 then
        -- n is quadratic residue (mod p) if Jacobi symbol (n / p) = 1
        return nil
    end

    if num.Words[1] % 4 == 3 then
        -- It is trivial to find the modular square root if p equiv 3 (mod 4).
        local modinv4 = module.modinv(4, p)
        return module.modexp(num, (p + 1) * modinv4, p)
    end

    -- Otherwise, we use the Tonelli-Shanks algorithm
    local totient_p = p - 1
    local S = totient_p:countrz()
    local Q = totient_p:shr(S)

    -- Finding a quadratic nonresidue mod p
    local z = BigInt:new(2)
    while module.nt.jacobi(z, p) ~= -1 do
        z += 1
    end

    local M = S
    local c = module.modexp(z, Q, p)
    local t = module.modexp(num, Q, p)
    local R = module.modexp(num, (Q + 1):shr(1), p)

    while true do
        if #t.Words == 1 then
            if t.Words[1] == 0 then
                return BigInt:new(0)
            end

            if t.Words[1] == 1 then
                return R
            end
        end

        local i = 1
        while i <= M do
            local exp = module.modexp(2, i, totient_p)

            if module.modexp(t, exp, p):eq(1) then
                break
            end

            i += 1
        end

        local E = module.modexp(2, M - i - 1, totient_p)
        local b = module.modexp(c, E, p)
        M = i
        c = (b * b) % p
        t = (t * c) % p
        R = (R * b) % p
    end
end

function module.sqrt (num: BigInt): BigInt
    --- Computes the integer square root i.e. floor(sqrt(n)) of a number.

    local x = module.internal.sqrt_optimization(num)

    local min_diff = num - (x ^ 2)
    local best_guess = x
    for i = 1, 5 do
        -- TODO: We can use the computation results to inform the next guess
        --       this can save significant time for extremely large numbers
        --       as the current implementation using naive algorithms is slow.
        local delta = BigInt:new(i)

        -- FIXME: luau-analyze really hates using symbols over the metamethods

        local xp = x + delta
        local xp2 = xp ^ 2

        if xp2:__eq(num) then
            return xp
        end

        local xm = x - delta
        local xm2 = xm ^ 2

        if xm2:__eq(num) then
            return xm
        end

        if xp2:__lt(num) then
            local diff = num - xp

            if (min_diff and diff:__lt(min_diff))
               or (not min_diff) then
                min_diff = diff
                best_guess = xp
            end
        end

        if xm2:__lt(num) then
            local diff = num - xm

            if (min_diff and diff:__lt(min_diff)) or (not min_diff) then
                diff = diff
                best_guess = xm
            end
        end
    end

    return best_guess
end

function module.nt.is_perfect_square (num: BigInt): boolean
    --- Checks whether a number is a perfect square.

    local x = module.internal.sqrt_optimization(num)
    for i = 1, 5 do
        -- TODO: See comment in equivalent loop in BigInt:sqrt()

        -- BigInt:sqrt() does additional checks at the end of the equivalent
        -- loop that will be unnecessarily slow when checking the squarness
        -- of very large numbers.
        local delta = BigInt:new(i)

        -- FIXME: LUAU-ANALYZE IS UNHAPPY WITH ==
        if ((x + delta) ^ 2):__eq(num) then
            return true
        end

        if ((x - delta) ^ 2):__eq(num) then
            return true
        end
    end

    return false
end

function module.nt.legendre (a: BigInt, p: BigInt): number
    -- Compute the Legendre symbol (a | p) where p prime.
    -- NOTE: We trust the user to know what they're doing using this function
    --       and gives a p prime and we avoid a costly primality check here.
    local remainder = a % p

    if remainder:eq(0) then
        return 0
    end

    local result = module.modexp(a, (p - 1) / 2, p)

    if result:eq(1) then
        return 1
    end

    return -1
end

function module.nt.jacobi2 (d: BigInt): number
    local m = d % 8

    if m:eq(1) or m:eq(7) then
        return 1
    end

    return -1
end

function module.nt.jacobi (a: BigInt, d: BigInt): number
    -- Computes the Jacobi symbol (a | d) where a positive and d negative,
    -- and the Kronecker symbol for all other a.
    a = a % d

    if a:eq(0) then
        return 0
    end

    local result = 1

    if a:eq(2) then
        return module.nt.jacobi2(d)
    end

    if a.Negative then
        if (d % 4):eq(3) then
            result = -1
        end

        a.Negative = false
    end

    while not a:eq(0) do
        if a.Words[1] % 2 == 0 then
            -- If a is even, then factor out powers of 2 such that a' := 2^r * a
            local r = a:countrz()
            a = a:shr(r)

            -- Jacobi (2 / d) is +/- 1, then (2 / d) ^ r for r even yields 1.
            -- Therefore, we only need to compute (2 / d) for r odd.
            if r % 2 > 0 then
                result *= module.nt.jacobi2(d)
            end
        end

        a, d = d, a
        if (a % 4):eq(3) and (d % 4):eq(3) then
            result = -result
        end

        a = a % d
    end

    if d:eq(1) then
        return result
    end

    return 0
end

function module.nt.internal.lucasmod_double (mod: BigInt,
                                             U: BigInt,
                                             V: BigInt,
                                             D: BigInt
                                            ): (BigInt, BigInt)
    local U_k = U

    U = (U_k * V) % mod
    local Vp = (V ^ 2) + D * (U_k ^ 2)
    if Vp.Words[1] % 2 == 1 then
        Vp += mod
    end
    V = Vp:shr(1) % mod

    return U, V
end

function module.nt.lucasmod (_num: Biggish,
                             _mod: Biggish,
                             _P: Biggish,
                             _D: Biggish
                            ): (BigInt, BigInt)
    -- This computes the nth number in the Lucas sequences
    -- U(P, Q) and V(P, Q) modulo m for use in Lucas probable prime test,
    -- where Q := (1 - D) / 4
    -- NOTE: We pass the discriminant D instead of Q to this function.
    local num: BigInt = BigInt:new(_num)
    local mod: BigInt = BigInt:new(_mod)
    local P: BigInt = BigInt:new(_P)
    local D: BigInt = BigInt:new(_D)

    local U = BigInt:new(0)
    local V = BigInt:new(2)

    for i = #num, 1, -1 do
        U, V = module.nt.internal.lucasmod_double(mod, U, V, D)

        if num:bit(i) == 1 then
            local U_2k = U

            local Up = P * U_2k + V
            if Up.Words[1] % 2 == 1 then
                Up += mod
            end
            U = Up:shr(1) % mod

            local Vp = D * U_2k + P * V
            if Vp.Words[1] % 2 == 1 then
                Vp += mod
            end
            V = Vp:shr(1) % mod
        end
    end

    return U, V
end

function module.nt.lucas (num: BigInt, P: BigInt, Q: BigInt): (BigInt, BigInt)
    -- Generate the nth numbers in the Lucas sequences U(P, Q) and V(P, Q)
    local D = P ^ 2 - 4 * Q
    local U = BigInt:new(0)
    local V = BigInt:new(2)

    for i = #num, 1, -1 do
        local U_k = U

        U = U * V
        V = ((V ^ 2) + D * (U_k ^ 2)):shr(1)

        if num:bit(i) == 1 then
            local U_2k = U

            U = (P * U_2k + V):shr(1)
            V = (D * U_2k + P * V):shr(1)
        end
    end

    return U, V
end

function module.prime.is_lucas_prime (num: BigInt): boolean
    --- Performs the regular Lucas probable prime test.
    --- We find a D such that the Jacobi (D / n) = -1, then compute the (n + 1)
    --- term in the Lucas sequence U(1, (1 - D) / 4). If U_{n + 1} is congruent
    --- to 0 modulo n, then it is a _probable_ prime.
    --- There are Lucas pseudoprimes (Lucas-Carmichael numbers) that pass this
    --- test, yet are composite.
    ---
    --- This method should _not_ be used, and instead one should consider using
    --- BigInt:is_strong_lucas_prime() instead.
    ---
    --- While all probablistic prime tests are prone to wrong answers,
    --- the regular Lucas test is susceptible to more false negatives than the
    --- strong Lucas test while requring about the same amount of computation.
    if num:eq(2) then
        return true
    end

    if (num.Words[1] % 2) == 0 then
        return false
    end

    local D = BigInt:new(5)

    local iters = 1
    while true do
        if iters == 12 then
            -- This method will loop infinitely if n is a perfect square.
            -- On average, it should take <= 3 steps to find a D such that the
            -- Jacobi (D / n) = -1, if it takes abnormally long then we should
            -- check if the number is a perfect square.

            -- The upperbound of 12 is arbitrarily chosen, and it should be
            -- exceptionally rare that we would call BigInt:is_perfect_square()
            -- which is good since it is considerably more expensive to perform
            -- Newton's than extraneous Jacobi checks.
            if module.nt.is_perfect_square(num) then
                return false
            end
        end
        local jacobi = module.nt.jacobi(D, num)

        if jacobi == -1 then
            break
        end

        if jacobi == 0 then
            -- If Jacobi symbol (D / num) = 0, then they share a common factor
            return false
        end

        local negative = not D.Negative
        D = D:abs() + 2
        D.Negative = negative

        iters += 1
    end

    local U, _ = module.nt.lucasmod(num + 1, num, 1, D)

    -- If U == 0, then it is a Lucas probable prime otherwise, it is composite.
    return U:eq(0)
end

function module.prime.is_strong_lucas_prime (num: BigInt): boolean
    --- Performs the strong Lucas primality test.
    if num:eq(2) or num:eq(3) or num:eq(5) then
        -- There's some weirdness when D = n, and this only happens when n = 5.
        -- So just bail and avoid doing the algorithm.
        return true
    end

    if (num.Words[1] % 2) == 0 then
        return false
    end

    local P = BigInt:new(3)
    local D = BigInt:new(5)

    local iters = 1
    while true do
        if iters == 12 then
            -- See comment for equivalent check in BigInt:lucas_prime()

            if module.nt.is_perfect_square(num) then
                return false
            end
        end
        D = (P ^ 2) - 4
        local jacobi = module.nt.jacobi(D, num)

        if jacobi == -1 then
            break
        end

        if jacobi == 0 then
            -- If Jacobi symbol (P / num) = 0, then they share a common factor
            return false
        end

        P += 1
    end

    local np1 = num + 1
    local s = np1:countrz()
    local d = np1:shr(s)
    local U_d, V_d = module.nt.lucasmod(d, num, P, D)

    if U_d:eq(0) then
        return true
    end

    -- If any V_{d + 2 ^ r} for r in {0, 1, 2, ..., s} is congruent to
    -- 0 mod n, then it is a strong Lucas probable prime
    if V_d:eq(0) then
        return true
    end

    local U_di, V_di = U_d, V_d
    for i = 1, s do
        U_di, V_di = module.nt.internal.lucasmod_double(num, U_di, V_di, D)

        if V_di:eq(0) then
            return true
        end
    end

    return false
end

function module.prime.is_miller_rabin_prime (num: BigInt, p_t: number?): boolean
    -- Checking whether a number is prime via Miller-Rabin

    if num:eq(2) then return true end
    if (num.Words[1] % 2) == 0 then return false end

    -- FIXME: We can optimize this.
    -- Sorenson and Webster (doi:10.1090/mcom/3134) show that testing with
    -- the first 14 primes as witnesses is adequate for < 82-bit integers
    -- Whereas bounding p_{k, t} according to FIPS 186-5 results in performing
    -- more tests for numbers < 82-bits than necessary.

    local np = num - 1
    local r = np:countrz()
    local d = np:shr(r)

    -- Always run an MR(2) test, it is the computationally cheapest test we
    -- can do, and serves as the basis for Baille-PSW test.
    local witnesses = { [2] = true }
    local k = #num

    if k <= 1024 then
        for i = 2, module.prime.internal.get_num_tests(k, p_t or 2^(-64)) do
            while true do
                local index = math.random(2, #module.prime.internal.primes)
                local prime = module.prime.internal.primes[index]

                if not witnesses[prime] then
                    witnesses[prime] = true
                    break
                end
            end
        end
    end

    for witness, _ in pairs(witnesses) do
        if module.prime.internal.miller_rabin(num, r, d, witness, np) then
            return false
        end
    end

    return true
end

function module.prime.is_prime (num: BigInt, p_t: number?): boolean
    --- Performs a conbination of Miller-Rabin with tolerance p_t and a strong
    --- Lucas probable prime test to determine whether a number is prime.
    ---
    --- Due to the usage of the strong Lucas probable prime test,
    --- we are able to use a smaller p_t for Miller-Rabin and thus fewer rounds
    --- of Miller-Rabin, accelerating the primality test.

    if num:eq(2) or num:eq(3) or num:eq(5) then
        -- Some weirdness happens when n=5 for the Lucas probable prime test.
        return true
    end

    if num.Words[1] % 2 == 0 then
        -- Is not 2 and is even.
        return false
    end

    if not module.prime.is_miller_rabin_prime(num, p_t or 2 ^ 16) then
        return false
    end

    return module.prime.is_strong_lucas_prime(num)
end

function module.nt.gcd (_a: Biggish, _b: Biggish): BigInt
    local a: BigInt = BigInt:new(_a)
    local b: BigInt = BigInt:new(_b)

    if a:eq(0) then
        return b:abs()
    end

    return module.nt.gcd(b:abs() % a:abs(), a:abs())
end

function module.nt.bezout_coefficients (_a: Biggish, _b: Biggish)
    -- Given a, b compute s, t such that s * a + t * b = (a, b)
    local a: BigInt = BigInt:new(_a)
    local b: BigInt = BigInt:new(_b)

    local pr = a:clone()
    local r = b:clone()

    local ps = BigInt:new(1)
    local s = BigInt:new(0)

    local pt = BigInt:new(0)
    local t = BigInt:new(1)

    while r:neq(0) do
        local q = pr / r

        pr, r = r, pr - (q * r)
        ps, s = s, ps - (q * s)
        pt, t = t, pt - (q * t)
    end

    return ps, pt
end

function module.modexp (_base: Biggish,
                        _exp: Biggish,
                        _mod: Biggish,
                        use_montgomery: boolean?
                       ): BigInt
    local base: BigInt = BigInt:new(_base)
    local exp: BigInt = BigInt:new(_exp)
    local mod: BigInt = BigInt:new(_mod)

    base = base % mod

    if base:eq(0) then
        return BigInt:new(0)
    end

    if base:eq(1) then
        return BigInt:new(1)
    end

    if base:eq(base - 1) then
        if exp.Words[1] % 2 == 1 then
            return base - 1
        end

        return BigInt:new(1)
    end

    -- TODO: Use exp := exp (mod totient(mod)) as an optimization
    --       However, computing totient(mod) needs fast factorization

    local accumulator = BigInt:new(1)
    if use_montgomery then
        error("modexp(use_montgomery=true) not implemented.")
    else
        for i = #exp, 1, -1 do
            accumulator = (accumulator * accumulator) % mod

            if exp:bit(i) == 1 then
                accumulator = (accumulator * base) % mod
            end
        end
    end

    return accumulator
end

function module.modinv (_num: Biggish, _mod: Biggish): BigInt
    local num: BigInt = BigInt:new(_num)
    local mod: BigInt = BigInt:new(_mod)

    local s, _ = module.nt.bezout_coefficients(num % mod, mod)
    return s % mod
end

return module
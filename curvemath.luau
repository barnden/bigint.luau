--!strict

local BigInt = require("./bigint")
local bmath = require("./bigmath")
local EC = require("./curves")
type BigInt = BigInt.BigInt
type Biggish = BigInt.Biggish
type Curve = EC.Curve
type Point = EC.Point

local module = {
    factorization = {
        internal = {}
    },
    signature = {
        internal = {},
        ecdsa = {
            internal = {},
        }
    }
}

function module.factorization.ecm (_n: Biggish): BigInt
    local n: BigInt = BigInt:new(_n)

    -- Generate a random Weierstrass curve over Z/nZ
    local a = bmath.random.random(#n + 1) % n
    local x0 = bmath.random.random(#n + 1) % n
    local y0 = bmath.random.random(#n + 1) % n
    local b = ((y0 ^ 2) - (x0 ^ 3) - a * x0) % n
    local curve: Curve = EC.Curve:new(n, a, b, EC.Types.Weierstrass)

    -- FIXME: Code too complex to use Point type here, for some reason
    local P: any = EC.Point:new(curve, x0, y0)
    local k = 1

    while true do
        local Q = k * P
        k += 1

        if Q.at_infinity then
            local d = bmath.nt.gcd(Q.x - P.x, n)

            if d:eq(1) or d:__eq(n) then
                return module.factorization.ecm(n)
            end

            return d
        end

        P = Q
    end
end

function module.signature.ecdsa.safe_sign (curve: Curve,
                                           G: Point,
                                           n: BigInt,
                                           d_A: BigInt,
                                           Q_A: BigInt,
                                           H_m: BigInt
                                          ): (BigInt, BigInt)
    -- Same as signature.ecdsa.sign, except we verify that n := ord(G) is prime
    -- and that n * G = O in E(Z/pZ). These are skipped in the other as checking
    -- these properties are nontrivial and incur high computational costs.

    if not bmath.prime.is_prime(n) then
        error("curvemath.signature.ecdsa.safe_sign() n := ord(G) is not prime.")
    end

    if not (G:__mul(n)).at_infinity then
        error(
            "curvemath.signature.ecdsa.safe_sign() n ~= ord(G) on "
            .. tostring(curve))
    end

    return module.signature.ecdsa.sign(curve, G, n, d_A, Q_A, H_m)
end

function module.signature.ecdsa.sign (curve: Curve,
                                      G: Point,
                                      n: BigInt,
                                      d_A: BigInt,
                                      Q_A: BigInt,
                                      H_m: BigInt
                                     ): (BigInt, BigInt)
    -- G: a base point on the curve such that ord(G) is prime
    -- n: ord(G) in curve
    -- d_A: private key
    -- Q_A: public key
    -- H_m: hashed message

    -- 1. Let z be L_n left-most bits of H_m, where L_n is the bit size of n
    --    H_m can be greater, but not longer, than n.
    local z: BigInt = H_m:clone()

    if #H_m > #n then
        z = z:shr(#H_m - #n)
    end

    local r: BigInt, s: BigInt
    repeat
        local k: BigInt, P: Point
        repeat
            -- 2. Select k ~= 0 in Z/nZ such that (k * G).x ~= 0
            -- NOTE: The PRNG is not a CSPRNG as required for this application.
            k = (bmath.random.random(#n) + 1) % n
            P = G * k
            r = P.x % n
        until k:neq(0)

        -- 3. Compute s, repeat if s is congruent to 0 mod n
        s = (bmath.modinv(k, n) * (z + r * d_A)) % n
    until s:neq(0)

    return r, s
end

function module.signature.ecdsa.verify (curve: Curve,
                                        G: Point,
                                        n: BigInt,
                                        Q_A: Point,
                                        H_m: BigInt,
                                        r: BigInt,
                                        s: BigInt
                                       ): boolean
    -- H_m: hashed message of which we wish to verify the signature
    -- r, s: the signature corresponding to H_m
    -- G: a base point on the curve such that ord(G) is prime
    -- n: ord(G) in curve
    -- Q_A: public key

    -- 1. Verify that Q_A is a valid curve point
    -- 1.1. Check that Q_A is not equal to the identity element
    if Q_A.at_infinity then
        return false
    end

    -- 1.2. Check that Q_A lies on the curve
    if not Q_A:on_curve() then
        return false
    end

    -- 1.3. Check that n is indeed ord(Q_A), i.e. Q_A * n = O in E[Z/pZ]
    if not (Q_A:__mul(n)).at_infinity then
        return false
    end

    -- 2. Verify that r, s are integers in [1, n - 1]
    if r:le(0) or s:le(0) then
        return false
    end

    if r:ge(n - 1) or s:ge(n - 1) then
        return false
    end

    -- 3. Calculate e = HASH(m); skipped as we assume H_m := HASH(m)
    -- 4. Let z be the L_n leftmost bits of H_m
    local z: BigInt = H_m:clone()

    if #H_m > #n then
        z = z:shr(#H_m - #n)
    end

    -- 5. Calculate u_1 = z * s ^ (-1) mod n and u_2 = r * s ^ (-1) mod n
    local s_inv = bmath.modinv(s, n)
    local u_1 = (z * s_inv) % n
    local u_2 = (r * s_inv) % n

    -- 6. Calculate the curve point P := (x_1, y_1) = u_1 * G + u_2 * Q_A
    local P: Point = (G:__mul(u_1) :: Point):__add(Q_A:__mul(u_2) :: Point)
    -- 6.1 If P = 0 then the signature is invalid
    if P.at_infinity then
        return false
    end

    -- 7. Accept if r equiv. x_1 (mod n).
    return r:__eq(P.x % n)
end

return module
# bigint.luau

An implementation of BigInt, and more, in Luau.

I suppose to teach the yutes.

## Usage
To create a `BigInt` call `BigInt:new(value: Biggish, radix: number?)`
- The type `Biggish` is the type union `BigInt | number | string`
- If `value` is a `BigInt`, it returns the _same_ object; _not_ a clone.
- If `value` is a `number` with magnitude $\leq 2^{53}$, otherwise `error()`.
- If `value` is a `string`, then the string is parsed in base-10 unless:
  - `radix` is specified; currently only `10` for decimal and `16` for hexadecimal are accepted, or
  - `value` is prefixed with `0x` denoting a hexadecimal string.

```lua
local a = BigInt:new(1337)
local b = BigInt:new("1234567890987654321")
local c = BigInt:new("d0229c8b57c5bd4cfd2693079d494c4f", 16)
local d = BigInt:new("0xd0229c8b57c5bd4cfd2693079d494c4f") -- same as above
```

`BigInt`s can be operated upon with the usual arithmetic operators, i.e.
```lua
local a: BigInt, b: BigInt

local      add = a + b
local subtract = a - b
local multiply = a * b
local    power = a ^ b
local   divide = a / b -- or a // b
local   modulo = a % b
```

Variables with type `Biggish` can also be used in these arighmetic operators,
```lua
local a: BigInt
local b = BigInt:new(5)

-- All of these yield equivalent results
local add = a + b
local add_number = a + 5
local add_string = a + "5"
```

The relational operators `<`, `<=`, `==`, `~=`, `=>`, `>` all work as expected between two `BigInt` objects, i.e.
```lua
local a = BigInt:new(29)
local b = BigInt:new(31)

print(a == b)         -- prints: false
print(a ~= b)         -- prints: true
print(a >= b)         -- prints: false
print((a + 2) == b)   -- prints: true
print(a == (b - "2")) -- prints: true
```

For comparing a `BigInt` with a `Biggish`, or a `Biggish` with a `Biggish`, the equivalent relational methods, defined below, must be used.

### Arithmetic and Relational Methods
- All arithmetic metamethods are implemented and accept `Biggish` arguments for the left- and right-hand side arguments.
  - `__div` and `__idiv` both perform integer division (floor).
  - `__mod` returns the modulo and not the remainder, i.e. `m % n` is always in $[0, n - 1]$
  - `__exp` does fast exponentation without using the Montgomery ladder method
    - the exponent is limited to a maximum of $2^{32} - 1$
    - negative exponents are invalid

- The `abs()` method returns a nonnegative valued clone of the number.

- Binary left (`shl`) and right (`shr`) bit-shifts are implemented.
  - Currently these do not accept `Biggish` left- and right-hand side arguments:
    - The left hand side is `BigInt` and right-hand is `BigInt | number`.
  - A bit-shift operation is is limited to $\leq 2 ^ {53}$ bits in either direction.

- The equality metamethods are also implemented, with equivalents which accept `Biggish` arguments for left- and right-hand side arguments.
  - `__eq` and `eq`
  - `__lt` and `lt`
  - `__le` and `le`
- Additional helper methods such as `neq`, `ge`, and `gt` are provided.


### Helper Methods
- The `clone()` method returns a deep-copy of a `BigInt`.
- The `__len` metamethod returns the number of binary digits.
- The `__tostring` metamethod returns a decimal representation of the `BigInt`.
- The `BigInt:hex(digits: number?)` method yields a hexadecimal string, without a `0x` prefix.
  - The `digits` parameter does left zero-padding on the hex string until reaching `digits` in length.

- The `bit(n: number)` method returns the `n`<sup>th</sup> binary digit, 1-indexed, where 1 corresponds to the least significant bit.
- The `countrz()` method returns the number of trailing zeros in the binary representation.
- The `value()` method returns a `number` literal in the range $(-2^{53}, 2^{53})$ with value equal to the `BigInt`, otherwise `nil`.

## BigMath
The `bigmath.luau` module defines some useful modular arithmetic, number theoretic, and random number generator methods.

### `BigMath`
The modular arithmetic operations are defined directly under the `BigMath` module

- `BigMath.modexp(base: Biggish, exp: Biggish, mod: Biggish): BigInt`
    - Computes `(base ^ exp) % mod`
    - Does _not_ use the Montgomery ladder method for exponentiation.
- `BigMath.modinv(num: Biggish, mod: Biggish): BigInt`
    - Computes the multiplicative inverse of $m \pmod{n}$ if it exists, i.e. `(num * modinv(num, mod)) % mod = 1`.
- `BigMath.modsqrt(num: Biggish, mod: Biggish): BigInt?`
    - Given a quadradic residue $q$ in $\mathbb{Z}/n\mathbb{Z}$, find $x$ such that  $x^2 \equiv q \pmod{n}$.
    - Returns `nil` if `num` is a quadradic nonresidue.
- `BigMath.sqrt(num: Biggish)`
    - Computes the integer square root of `num`.

### `BigMath.random`
- `BigMath.random.random(bits: number): BigInt`
    - Generates a random `BigInt` with `<= bits` binary digits.
    - Relies on the underlying Luau `random.random(a, b)` PRNG.
- `BigMath.random.random_prime(bits: number): BigInt`
    - Generates a random odd prime with `<= bits` binary digits.

### `BigMath.nt`
- `BigMath.nt.lucas(num: BigInt, P: BigInt, Q: BigInt)`
    - Generates the n<sup>th</sup> values in the Lucas sequences $U(P, Q)$ and $V(P, Q)$.
    - Note that each step does not have constant work.
- `BigMath.nt.lucasmod(num: Biggish, mod: Biggish, P: Biggish, D: Biggish)`
    - Generates the n<sup>th</sup> values in the Lucas sequences $U(P, Q)$ and $V(P, Q)$ modulo $m$, where $Q \triangleq (1 - D) / 4$.
    - Note that each step does not have constant work.
- `BigMath.nt.jacobi(a: BigInt, d: BigInt): number`
    - Computes the Jacobi symbol $(a / d)$
- `BigMath.nt.jacobi2(d: BigInt): number`
    - Computes the Jacobi symbol $(2 / d)$
- `BigMath.nt.legendre(a: BigInt, p: BigInt): number`
    - Computes the Legendre symbol $(a / p)$ where $p$ prime
    - The function assumes the user knows what they're doing and does not verify $p$ prime.
- `BigMath.nt.is_perfect_square(num: BigInt): boolean`
    - Checks whether a number is a perfect square, equivalent to `BigMath.sqrt(num) == num`.
- `BigMath.nt.gcd(a: Biggish, b: Biggish)`
    - Computes the greatest common divisor between $a$ and $b$.
- `BigMath.nt.bezout_coefficients(a: Biggish, b: Biggish): (BigInt, BigInt)`
    - Finds the Bezout coefficients $s, t$ such that $sa + bt = \gcd(a, b)$

### `BigMath.prime`
- `BigMath.prime.is_lucas_prime(num: BigInt): boolean`
    - Checks if a number is prime, using the Lucas probable prime test
- `BigMath.prime.is_strong_lucas_prime(num: BigInt): boolean`
    - Checks if a number is prime, using the strong Lucas probable prime test
- `BigMath.prime.is_miller_rabin_prime(num: BigInt, p_t: number?): boolean`
    - Checks if a number is prime, using Miller-Rabin
    - `p_t` is the threshold probability (default: $2^{-64}$) to determine the
      number of witnesses we should test against, the number of tests is determined
      by the formula given in FIPS 186-5.
- `BigMath.prime.is_prime(num: BigInt, p_t: number?): boolean`
    - Checks if a number is prime, effectively performing a Baille-PSW test.
    - `p_t` is the threshold probability (see above and FIPS 186-5) defaulting to a $2 ^ {-16}$ as combining Miller-Rabin with a strong Lucas probable prime test means that it extremely unlikely a number is strongly pseudoprime to both tests.

## Elliptic Curves
A basic elliptic curve module is provided in `curves.luau`.

### Usage
To define a curve, call `Curve:new(field: Biggish, a: Biggish, b: Biggish, type: CurveType?)`
- Where `field` is a prime $p$ and $a, b \in \mathbb{Z}/p\mathbb{Z}$.
- A curve can either be in Weierstrass or Montgomery form by specifying `type` as `CurveTypes.Weierstrass` and `CurveTypes.Montgomery`, respectively.
    - By default, `type` uses `CurveTypes.Weierstrass`
    - In the module, `CurveTypes` is exported as `Types`
- A curve given in Montgomery form is converted into the isomorphic Weierstrass curve internally.
    - Some algorithms assume Montgomery form for more efficient group operations, but we lose this due to the internal Weierstrass representation.

To evaluate the curve, invoke the `__call` metamethod, i.e. let `curve = my_curve:new(...)` then `point: Point = my_curve(x: Biggish)`.
  - The `x` value is given in the form of the original curve, i.e. if we created a Montgomery curve we need not convert the `x` into Weierstrass form.
  - Internally, the resulting `Point` object uses Weierstrass form coordinates, however, `__tostring` will print out the coordinates in the curve's original form.

The arithmetic operators `+` and `*` are defined as follows:
- `__add(a: Point, b: Point): Point` is defined to be the elliptic curve group operation $a + b$ between two points.
- `__mul(a: Point | Biggish, b: Point | Biggish): Point` is defined as the repeated addition of a point onto itself.
    - Again, this does not use the Montgomery ladder method.

### Standard Curves
A `StandardCurves` table is added for convenience defining a couple of commonly used elliptic curves.

Currently included curves are:
- `secp256k1`
- `secp256r1` or `NIST_P256`
- `Curve25519`

Each curve has the following properties
- `Curve` the curve object which can be evaluated to obtain a point,
- `Type` the form of the curve,
- `Generator` is the generator point,
- `n` is the order of the generator point,
- `p` the prime field,
- `a` and `b` are parameters defining the curve,
- `g` is the $x$-coordinate used in the generator

### Example
We can define the popular Curve25519 as follows:
```lua
local EC = require("./curves.luau")
local p25519 = (BigInt:new(2) ^ 255) - 19
local Curve25519 = EC.Curve:new(p25519, 486662, 1, EC.CurveTypes.Montgomery)

-- Generate the typical base point at x=9 (in Montgomery coordinates)
local G = Curve25519(9)
-- G := Point(x=9, y=14781619447589544791020593568409986887264606134616475288964881837755586237401)
```

## CurveMath
The CurveMath module provided by `curvemath.luau` provides some functions that operate on elliptic curves.

- `CurveMath.shamir(u: Biggish, G: Point, v: Biggish, Q: Point)`
  - Computes $u \cdot G + v \cdot Q$ where $u, v \in \mathbb{Z}$ and $G, Q \in E[\mathbb{F}_p]$ using Shamir's (aka Straus') trick

### `CurveMath.factorization`
  - `CurveMath.factorization.ecm(n: Biggish): BigInt`
    - Generates a random Weierstrass curve $E[\mathbb{Z}/n\mathbb{Z}]$ and performs Lenstra ellitpic curve factorization in order to find one random factor of $n$.
    - Assumption is that the user knows what they're doing and does not pass an $n$ that is already prime.
    - Due to randomness, run-time is highly variable and dependent on the size of the prime factors.
### `CurveMath.signature.ecdsa`
Compute and verify ellitpic curve digital signatures, make sure to bring your own hash function!

  - `CurveMath.signature.ecdsa.safe_sign(curve: Curve, G: Point, n: BigInt, d_A: BigInt, Q_A: BigInt, H_m: BigInt): (BigInt, BigInt)`
    - Same as `CurveMath.signature.ecdsa.sign` except we check $n \stackrel{?}{=} \mathrm{ord}(G)$ and that $n \cdot G \in E[\mathbb{Z}/p\mathbb{Z}]$
        - These checks are nontrivial and incur high compute costs; as such, they are skipped normally as we assume the caller to know what they're doing.
  - `CurveMath.signature.ecdsa.sign(curve: Curve, G: Point, n: BigInt, d_A: BigInt, Q_A: BigInt, H_m: BigInt): (BigInt, BigInt)`
    - Computes the ECDSA signature $(r, s)$ using the elliptic curve $E[\mathbb{Z}/p\mathbb{Z}]$ where
      - `G` is a base point on $E[\mathbb{Z}/p\mathbb{Z}]$ such that $\mathrm{ord}(G)$ is prime,
      - `n` is the order of the base point $G$ (see above),
      - `d_A` is the private key which is a randomly chosen element of $E[\mathbb{Z}/p\mathbb{Z}]$,
      - `Q_A` is the public key defined as the point $G \cdot d_A$,
      - `H_m` is the hashed message
  - `CurveMath.signature.ecdsa.ASN1DER(r: BigInt, s: BigInt): string`
    - Converts the public key from an ordered pair $(r, s)$ into a ASN.1/DER encoding
  - `CurveMath.signature.ecdsa.verify(curve: Curve, G: Point, n: BigInt, Q_A: Point, H_m: BigInt, r: BigInt, s: BigInt): boolean`
    - Determines whether the provided $(r, s)$ constitute a valid ECDSA signature for the hashed message `H_m` on the given curve.
      - `G` is a base point on $E[\mathbb{Z}/p\mathbb{Z}]$ such that $\mathrm{ord}(G)$ is prime,
      - `n` is the order of the base point $G$ (see above),
      - `Q_A` is the public key, 
      - `H_m` is the hashed message,
      - `r` and `s` is the ECDSA signature

### Example ECDSA
```lua
local BigInt = require("./bigint")
local bmath = require("./bigmath")
local EC = require("./curves")
local cmath = require("./curvemath")

math.randomseed(0)

-- Using secp256k1 curve defined in curves.luau
local curve = EC.StandardCurves.secp256k1

local secp256k1 = curve.Curve
local G = curve.Generator -- The generator point
local ord_G = curve.n -- Order of the generator point on the curve

-- Generate a random private key from [1, n - 1]
local d_A = bmath.random.random_modulo(curve.n, true)
print("ECDSA privkey ........ (hex): " .. d_A:hex())
-- prints: ECDSA privkey ........ (hex): 3b6eb58b2ebad58fa4c17181d755886f3114dba032dbb5504496bd00bcf68726

local length = math.ceil(#ord_G / 4) -- used for padding hex strings in output

-- Generate public key using private key and generator point
local Q_A = G * d_A
print("ECDSA pubkey .. (ANSI X9.63): 04" .. Q_A.x:hex(length) .. Q_A.y:hex(length))
-- prints: ECDSA pubkey .. (ANSI X9.63): 04282ee9e14eecc9723ea6fe9eedbcb8c770fb8b38b0709b0cde6040d0a5f6be320f716dfb474630bfc52a4da2bbe0cd058fcea3846196ceed6bf3b8bdaf1b6afc

-- The hashed message, in this case H_m := SHA256("Hello, World!") using utf-8
local H_m = BigInt:new("0xdffd6021bb2bd5b0af676290809ec3a53191dd81c7f70a4b28688a362182986f")

-- Generate signature for H_m using secp256k1
local r, s = cmath.signature.ecdsa.sign(secp256k1, G, ord_G, d_A, Q_A, H_m)

print("ECDSA Signature ..... (r||s): " .. r:hex(length) .. s:hex(length))
-- prints: ECDSA Signature ..... (r||s): d044c8fae940af1d99370b73e7960cf46d6c33c4934aaaa15069ddd3698a909f04179d6230a1359207541916817b4afa0d3ce28523196e59ce823a88e108bf2b
print("ECDSA Signature  (ASN.1/DER): " .. cmath.signature.ecdsa.ASN1DER(r, s))
-- prints: ECDSA Signature  (ASN.1/DER): 3045022100d044c8fae940af1d99370b73e7960cf46d6c33c4934aaaa15069ddd3698a909f022004179d6230a1359207541916817b4afa0d3ce28523196e59ce823a88e108bf2b

-- Check that the signature we just created is valid
local verify = cmath.signature.ecdsa.verify(secp256k1, G, ord_G, Q_A, H_m, r, s)
print("ECDSA Signature    verified?: " .. tostring(verify))
-- prints: ECDSA Signature    verified?: true
```
## Implementation Details
- `BigInt`s are represented as an array of base $2^{16}$ digits with a `boolean` defining whether a number is negative or nonegative.
    - This is due to the fact that Luau does not support 64-bit integers, so we rely on the `bit32` library to perform bitwise operations.
- Multiplication is performed using standard, or grade-school, multiplication up until numbers around $2^{768}$
  - After this limit, Karatsuba's multiplication method is used for the algorithmic improvement.
- Division is performed using standard division, as described in Algorithm D from _The Art of Computer Programming, Vol. 2_ by Donald Knuth.
- Exponentiation is performed using the fast powering method, not using the Montgomery ladder.


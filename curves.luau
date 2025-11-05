--!strict

local BigInt = require("./bigint")
local bmath = require("./bigmath")
type BigInt = BigInt.BigInt
type Biggish = BigInt.Biggish

type CurveTypeData = {
    Unique: { any },
    Name: string
}
local CurveType = {}
CurveType.__index = CurveType
function CurveType:new(name): CurveType
    local instance = { Unique = {}, Name = name }
    setmetatable(instance, self)

    return instance
end

function CurveType.__tostring (Type: CurveType): string
    return Type.Name
end

export type CurveType = typeof(setmetatable({} :: CurveTypeData, CurveType))

local CurveTypes = {
    Montgomery = CurveType:new("Montgomery"),
    Weierstrass = CurveType:new("Weierstrass"),
    -- TwistedEdwards = CurveType:new(),
}

type CurveData = {
    Type: CurveType,
    Field: BigInt,
    a: BigInt,
    b: BigInt,
    x_to_weierstrass: (BigInt) -> BigInt,
    y_to_weierstrass: (BigInt) -> BigInt,
    x_from_weierstrass: (BigInt) -> BigInt,
    y_from_weierstrass: (BigInt) -> BigInt,
}
local Curve = {}
Curve.__index = Curve
export type Curve = typeof(setmetatable({} :: CurveData, Curve))

local Point = {}
Point.__index = Point
type PointData = {
    Curve: Curve,
    x: BigInt,
    y: BigInt,
    at_infinity: boolean,
}
export type Point = typeof(setmetatable({} :: PointData, Point))

function __characterize (curve: Curve): (BigInt, BigInt?)
    local c4 = (-48 * curve.a) % curve.Field
    local c6 = (-864 * curve.b) % curve.Field

    local modinv1728 = bmath.modinv(1728, curve.Field)
    local DeltaE = ((c4 ^ 3 - c6 ^ 2) * modinv1728) % curve.Field

    if DeltaE:eq(0) then
        -- Curve is singular
        return DeltaE, nil
    end

    local invDeltaE = bmath.modinv(DeltaE, curve.Field)
    local j_invariant = ((c4 ^ 3) * invDeltaE) % curve.Field

    return DeltaE, j_invariant
end

function Curve:new (_field: Biggish,
                    _a: Biggish,
                    _b: Biggish,
                    _Type: CurveType?
                   ): Curve
    --- Construct a new Weierstrass elliptic curve.
    --- Other curve types are converted into Weierstrass form.
    ---
    --- Points generated on the curve printed out via tostring() convert the
    --- coordinates back to their original form, rather than in Weierstrass form
    --- The actual point data, i.e. members x and y, are in Weierstrass form
    local field: BigInt = BigInt:new(_field)
    local a: BigInt = BigInt:new(_a)
    local b: BigInt = BigInt:new(_b)

    local Type: CurveType = _Type or CurveTypes.Weierstrass

    local instance = setmetatable({} :: CurveData, Curve)

    instance.Type = Type
    instance.Field = field
    instance.a = a
    instance.b = b

    instance.x_to_weierstrass = function (x: BigInt): BigInt
        return x
    end
    instance.y_to_weierstrass = function (y: BigInt): BigInt
        return y
    end
    instance.x_from_weierstrass = function (x: BigInt): BigInt
        return x
    end
    instance.y_from_weierstrass = function (y: BigInt): BigInt
        return y
    end

    if Type == CurveTypes.Montgomery then
        if b:eq(1) then
            local dx = (a * bmath.modinv(3, field)) % field

            instance.x_to_weierstrass = function (x: BigInt): BigInt
                return (x + dx) % field
            end

            instance.y_to_weierstrass = function (y: BigInt): BigInt
                return y % field
            end

            instance.x_from_weierstrass = function (x: BigInt): BigInt
                return (x - dx) % field
            end

            instance.y_from_weierstrass = function (y: BigInt): BigInt
                return y % field
            end
        else
            local modinvB = bmath.modinv(b, field)
            local dx = (a * bmath.modinv(3 * b, field)) % field

            instance.x_to_weierstrass = function (x: BigInt): BigInt
                return (x * modinvB + dx) % field
            end

            instance.y_to_weierstrass = function (y: BigInt): BigInt
                return (y * modinvB) % field
            end

            instance.x_from_weierstrass = function (x: BigInt): BigInt
                return ((x - dx) * b) % field
            end

            instance.y_from_weierstrass = function (y: BigInt): BigInt
                return (y * b) % field
            end
        end

        instance.a = (3 - (a ^ 2)) * bmath.modinv(3 * b ^ 2, field)
        instance.b = (2 * (a ^ 3) - (9 * a)) * bmath.modinv(27 * (b ^ 3), field)

        instance.a %= field
        instance.b %= field
    end

    -- FIXME: adding discriminant to CurveData causes luau-analyze errors
    local _discriminant
    _discriminant, instance.j_invariant = __characterize(instance :: Curve)

    return table.freeze(instance)
end

function Curve.__tostring (curve: Curve): string
    return "Curve("
           .. "Type=" .. tostring(curve.Type) .. ", "
           .. "Field=" .. tostring(curve.Field) .. ", "
           .. "a=" .. tostring(curve.a) .. ", "
           .. "b=" .. tostring(curve.b)
           ..")"
end

function Curve.__call (curve: Curve, _x: Biggish): Point?
    local x: BigInt = BigInt:new(_x)

    local t: BigInt = curve.x_to_weierstrass(x % curve.Field)
    -- FIXME: luau-analyze does not recognize correct type with '%'
    local y2: BigInt = ((t ^ 3) + curve.a * t + curve.b):__mod(curve.Field)
    local y: BigInt? = bmath.modsqrt(y2, curve.Field)

    if not y then
        return nil
    end

    local v: BigInt = curve.y_to_weierstrass(y)
    return Point:new(curve, t, v)
end

function Point:new (curve: Curve,
                    x: BigInt,
                    y: BigInt,
                    at_infinity: boolean?
                   ): Point
    -- Defines a homogeneous coordinate along a curve, where w denotes whether
    -- a point lies at infinity, or on the curve.
    --
    -- Regardless of curve type, the coordinates of the point are defined w.r.t
    -- the Weierstrass form.
    --
    -- When printed via Point.__tostring(), the coordinates are displayed in
    -- the form for the curve type, and not in Weierstrass form.
    local instance = setmetatable({} :: PointData, Point)

    instance.Curve = curve
    instance.x = x
    instance.y = y
    instance.at_infinity = at_infinity == true

    return instance
end

function Point:copy (point: Point)
    self.Curve = point.Curve
    self.x = point.x
    self.y = point.y
    self.at_infinity = point.at_infinity
end

function Point.clone (point: Point): any
    -- FIXME: Code is too complex to analyze if we don't explicitly pass Point
    return Point.new(Point, point.Curve, point.x, point.y, point.at_infinity)
end

function Point.infinity (point: Point): Point
    local result = point:clone()

    result.x = BigInt:new(0)
    result.y = BigInt:new(0)
    result.at_infinity = true

    return result
end

function Point.__eq (lhs: Point, rhs: Point): boolean
    return (lhs.Curve == rhs.Curve) and (lhs.x == rhs.x) and (lhs.y == rhs.y)
end

function Point.__unm (point: Point): Point
    local result = point:clone()

    result.y = result.Curve.Field - result.y

    return result
end

function Point.__add (lhs: Point, rhs: Point): Point
    if (lhs.Curve ~= rhs.Curve) then
        -- FIXME: Needs more robust error handling
        return lhs:infinity()
    end

    if lhs.at_infinity and rhs.at_infinity then
        return lhs:infinity()
    end

    if lhs.at_infinity then
        return rhs:clone()
    end

    if rhs.at_infinity then
        return lhs:clone()
    end

    local curve = lhs.Curve
    local xp = lhs.x
    local yp = lhs.y

    local xq = rhs.x
    local yq = rhs.y

    if (yp ~= yq) and (xp == xq) then
        -- Points lie on a vertical line, ergo ad infinitum
        -- NOTE: Could we also treat this as a degenerate third point which
        -- overlaps one of the existing two points? Probably won't be useful
        -- in this application, however.

        return lhs:infinity()
    end

    local lambda
    -- FIXME: luau-analyze code too complex if using == instead of __eq()
    if lhs:__eq(rhs) then -- Specialize for point doubling
        lambda = (3 * (xp ^ 2) + curve.a) * bmath.modinv(2 * yp, curve.Field)
    else
        lambda = (yq - yp) * bmath.modinv(xq - xp, curve.Field)
    end
    lambda %= curve.Field

    local xr = ((lambda ^ 2) - xp - xq) % curve.Field
    local yr = ((lambda * (xp - xr) - yp)) % curve.Field

    local expected_y2 = (xr ^ 3 + curve.a * xr + curve.b) % curve.Field
    local computed_y2 = (yr ^ 2) % curve.Field

    if not computed_y2:__eq(expected_y2) then
        return lhs:infinity()
    end

    return Point:new(curve, xr, yr)
end

function Point.__sub (lhs: Point, rhs: Point): Point
    -- TODO: Handle when lhs and rhs are not on the same curve
    if rhs.at_infinity then
        return lhs
    end

    -- FIXME: luau-analyze code too complex if using + instead of __add()
    return lhs:__add(-rhs)
end

function Point.__mul (lhs: Point | Biggish, rhs: Point | Biggish): Point
    local point: Point
    local factor: BigInt

    if type(lhs) == "table" and getmetatable(lhs) == Point then
        point = lhs :: Point
        factor = BigInt:new(rhs :: Biggish)
    elseif type(rhs) == "table" and getmetatable(rhs) == Point then
        point = rhs :: Point
        factor = BigInt:new(lhs :: Biggish)
    else
        error("TypeError")
    end

    if factor:eq(0) then
        return point:infinity()
    end

    if factor:lt(0) then
        point = -point
    end

    factor = factor:abs()
    if factor:eq(1) then
        return point
    end

    -- Point at infinity acts as our identity element under addition
    local result: any = point:infinity()

    -- We should clamp factor to a u32
    for i = #factor, 1, -1 do
        -- FIXME: luau-analyze code too complex to analyze if using +=
        result = result:__add(result)

        if factor:bit(i) == 1  then
            -- FIXME: luau-analyze code too complex to analyze if using +=
            result = result:__add(point)
        end
    end

    return result
end

function Point.on_curve (point: Point): boolean
    if point.at_infinity then return true end

    local curve = point.Curve

    local x = curve.x_from_weierstrass(point.x)
    -- FIXME: luau-analyze code too complex to analyze if using curve(x)
    local P = curve:__call(x)

    if not P then
        return false
    end

    local computed = bmath.modexp(point.y, 2, curve.Field)
    local expected = bmath.modexp(P.y, 2, curve.Field)
    return computed:__eq(expected)
end

function Point.__tostring (point: Point): string
    if point.at_infinity then
        return "Point(at_infinity=true)"
    end

    return "Point("
         .. "x=" .. tostring(point.Curve.x_from_weierstrass(point.x)) .. ", "
         .. "y=" .. tostring(point.Curve.y_from_weierstrass(point.y))
         .. ")"
end

return {
    ["Curve"] = Curve,
    ["Point"] = Point,
    ["Types"] = CurveTypes
}
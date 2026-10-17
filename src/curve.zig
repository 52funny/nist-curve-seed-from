const std = @import("std");
const BigInt = std.math.big.int.Managed;

/// Solves c*b^2 = -27 (mod p) for the curve y^2 = x^3 - 3*x + b.
/// p must be a prime greater than 3; primality is a caller precondition.
/// Returns a root in [0, p); the other root is p-b. Caller must deinit it.
/// For p = 3 (mod 4), selects (-27/c)^((p+1)/4), as in the reference.
/// Rejects zero c, singular curves, and coefficients without a square root.
pub fn curveBFromC(allocator: std.mem.Allocator, c: *const BigInt, p: *const BigInt) !BigInt {
    if (p.toConst().orderAgainstScalar(3) != .gt or !p.isOdd())
        return error.InvalidPrime;

    var field = try PrimeField.init(allocator, p);
    defer field.deinit();
    var normalized = try BigInt.init(allocator);
    defer normalized.deinit();
    try field.reduce(&normalized, c);
    if (normalized.eqlZero()) return error.InvalidCoefficient;

    var check = try BigInt.init(allocator);
    defer check.deinit();
    try check.shiftLeft(&normalized, 2);
    try check.addScalar(&check, 27);
    try field.reduce(&check, &check);
    if (check.eqlZero()) return error.SingularCurve;

    var inverse = try field.inverse(&normalized);
    defer inverse.deinit();

    var minus27 = try BigInt.initSet(allocator, -27);
    defer minus27.deinit();
    var squared_b = try BigInt.init(allocator);
    defer squared_b.deinit();
    try field.mul(&squared_b, &minus27, &inverse);

    var b = try field.sqrt(&squared_b);
    errdefer b.deinit();
    try field.mul(&check, &b, &b);
    if (!check.eql(squared_b)) return error.NoSquareRoot;
    return b;
}

/// Computes c = -27/(b*b) mod p once when searching for a known coefficient b.
/// p must be prime and greater than 3. Caller must deinit the result.
/// Both b and p-b produce the same c; root selection must be checked separately.
pub fn curveCFromB(allocator: std.mem.Allocator, b: *const BigInt, p: *const BigInt) !BigInt {
    if (p.toConst().orderAgainstScalar(3) != .gt or !p.isOdd())
        return error.InvalidPrime;

    var field = try PrimeField.init(allocator, p);
    defer field.deinit();
    var squared_b = try BigInt.init(allocator);
    defer squared_b.deinit();
    try field.mul(&squared_b, b, b);
    if (squared_b.eqlZero()) return error.InvalidCoefficient;
    // With a = -3, the curve is singular exactly when b*b = 4 mod p.
    if (squared_b.toConst().orderAgainstScalar(4) == .eq)
        return error.SingularCurve;

    var inverse = try field.inverse(&squared_b);
    defer inverse.deinit();
    var minus27 = try BigInt.initSet(allocator, -27);
    defer minus27.deinit();
    var c = try BigInt.init(allocator);
    errdefer c.deinit();
    try field.mul(&c, &minus27, &inverse);
    return c;
}

// All methods allow the output to alias their inputs, except the modulus.
const PrimeField = struct {
    allocator: std.mem.Allocator,
    p: *const BigInt,
    quotient: BigInt,
    product: BigInt,

    fn init(allocator: std.mem.Allocator, p: *const BigInt) !PrimeField {
        var quotient = try BigInt.init(allocator);
        errdefer quotient.deinit();
        return .{
            .allocator = allocator,
            .p = p,
            .quotient = quotient,
            .product = try BigInt.init(allocator),
        };
    }

    fn deinit(self: *PrimeField) void {
        self.quotient.deinit();
        self.product.deinit();
    }

    fn reduce(self: *PrimeField, out: *BigInt, value: *const BigInt) !void {
        try BigInt.divFloor(&self.quotient, out, value, self.p);
    }

    fn mul(self: *PrimeField, out: *BigInt, a: *const BigInt, b: *const BigInt) !void {
        try self.product.mul(a, b);
        try self.reduce(out, &self.product);
    }

    fn inverse(self: *PrimeField, value: *const BigInt) !BigInt {
        // Fermat's little theorem: value^(-1) = value^(p-2) mod p.
        var exponent = try self.p.cloneWithDifferentAllocator(self.allocator);
        defer exponent.deinit();
        try exponent.addScalar(&exponent, -2);
        var result = try self.pow(value, &exponent);
        errdefer result.deinit();
        var check = try BigInt.init(self.allocator);
        defer check.deinit();
        try self.mul(&check, value, &result);
        if (check.toConst().orderAgainstScalar(1) != .eq)
            return error.InvalidCoefficient;
        return result;
    }

    fn pow(self: *PrimeField, base: *const BigInt, exponent: *const BigInt) !BigInt {
        var result = try BigInt.initSet(self.allocator, 1);
        errdefer result.deinit();
        var factor = try BigInt.init(self.allocator);
        defer factor.deinit();
        try self.reduce(&factor, base);
        var e = try exponent.cloneWithDifferentAllocator(self.allocator);
        defer e.deinit();

        while (!e.eqlZero()) {
            if (e.isOdd()) try self.mul(&result, &result, &factor);
            try e.shiftRight(&e, 1);
            if (!e.eqlZero()) try self.mul(&factor, &factor, &factor);
        }
        return result;
    }

    fn sqrt(self: *PrimeField, n: *const BigInt) !BigInt {
        var q = try self.p.cloneWithDifferentAllocator(self.allocator);
        defer q.deinit();
        try q.addScalar(&q, -1);
        var s: usize = 0;
        while (!q.isOdd()) : (s += 1) try q.shiftRight(&q, 1);

        var exponent = try BigInt.init(self.allocator);
        defer exponent.deinit();
        if (s == 1) {
            // Fast path for P-192, P-256, P-384, and P-521.
            try exponent.addScalar(self.p, 1);
            try exponent.shiftRight(&exponent, 2);
            return self.pow(n, &exponent);
        }

        // Tonelli-Shanks handles primes p = 1 (mod 4), including P-224.
        try exponent.addScalar(self.p, -1);
        try exponent.shiftRight(&exponent, 1);
        var symbol = try self.pow(n, &exponent);
        defer symbol.deinit();
        if (symbol.toConst().orderAgainstScalar(1) != .eq)
            return error.NoSquareRoot;

        var minus_one = try self.p.cloneWithDifferentAllocator(self.allocator);
        defer minus_one.deinit();
        try minus_one.addScalar(&minus_one, -1);
        var z = try BigInt.initSet(self.allocator, 2);
        defer z.deinit();
        while (z.order(self.p.*) == .lt) {
            var z_symbol = try self.pow(&z, &exponent);
            defer z_symbol.deinit();
            if (z_symbol.eql(minus_one)) break;
            try z.addScalar(&z, 1);
        } else return error.NoSquareRoot;

        var factor = try self.pow(&z, &q);
        defer factor.deinit();
        var t = try self.pow(n, &q);
        defer t.deinit();
        try exponent.addScalar(&q, 1);
        try exponent.shiftRight(&exponent, 1);
        var root = try self.pow(n, &exponent);
        errdefer root.deinit();
        var m = s;

        while (t.toConst().orderAgainstScalar(1) != .eq) {
            var probe = try t.clone();
            defer probe.deinit();
            var i: usize = 0;
            while (i < m and probe.toConst().orderAgainstScalar(1) != .eq) : (i += 1)
                try self.mul(&probe, &probe, &probe);
            if (i == m) return error.NoSquareRoot;

            var step = try factor.clone();
            defer step.deinit();
            for (0..m - i - 1) |_| try self.mul(&step, &step, &step);
            try self.mul(&root, &root, &step);
            try self.mul(&factor, &step, &step);
            try self.mul(&t, &t, &factor);
            m = i;
        }
        return root;
    }
};

test "curveBFromC agrees with exhaustive search over small prime fields" {
    for ([_]u32{ 5, 7, 11, 13, 17, 19, 23, 29, 41, 97 }) |prime| {
        var p = try BigInt.initSet(std.testing.allocator, prime);
        defer p.deinit();
        var c = try BigInt.init(std.testing.allocator);
        defer c.deinit();
        for (0..prime + 2) |value| {
            try c.set(value);
            const reduced: u32 = @intCast(value % prime);
            if (reduced == 0) {
                try std.testing.expectError(error.InvalidCoefficient, curveBFromC(std.testing.allocator, &c, &p));
                continue;
            }
            if ((4 * reduced + 27) % prime == 0) {
                try std.testing.expectError(error.SingularCurve, curveBFromC(std.testing.allocator, &c, &p));
                continue;
            }
            var has_root = false;
            for (0..prime) |candidate| {
                if ((reduced * candidate * candidate + 27) % prime == 0) {
                    has_root = true;
                    break;
                }
            }
            if (!has_root) {
                try std.testing.expectError(error.NoSquareRoot, curveBFromC(std.testing.allocator, &c, &p));
                continue;
            }
            var b = try curveBFromC(std.testing.allocator, &c, &p);
            defer b.deinit();
            const actual = try b.toInt(u32);
            try std.testing.expect(actual < prime);
            try std.testing.expectEqual(@as(u32, 0), (reduced * actual * actual + 27) % prime);
            var recovered_c = try curveCFromB(std.testing.allocator, &b, &p);
            defer recovered_c.deinit();
            try std.testing.expectEqual(reduced, try recovered_c.toInt(u32));
        }
    }
}

test "curveBFromC normalizes negative coefficients and rejects invalid moduli" {
    var p = try BigInt.initSet(std.testing.allocator, 41);
    defer p.deinit();
    var c = try BigInt.initSet(std.testing.allocator, -38);
    defer c.deinit();
    var b = try curveBFromC(std.testing.allocator, &c, &p);
    defer b.deinit();
    const actual = try b.toInt(u32);
    try std.testing.expect(actual == 14 or actual == 27);

    for ([_]i32{ -7, 0, 1, 2, 3, 4, 10 }) |invalid| {
        try p.set(invalid);
        try std.testing.expectError(error.InvalidPrime, curveBFromC(std.testing.allocator, &c, &p));
    }
}

fn checkAllocationFailures(allocator: std.mem.Allocator) !void {
    var p = try BigInt.initSet(allocator, 41);
    defer p.deinit();
    var c = try BigInt.initSet(allocator, 3);
    defer c.deinit();
    var b = try curveBFromC(allocator, &c, &p);
    defer b.deinit();
    const actual = try b.toInt(u32);
    try std.testing.expect(actual == 14 or actual == 27);
    var recovered_c = try curveCFromB(allocator, &b, &p);
    defer recovered_c.deinit();
    try std.testing.expect(c.eql(recovered_c));
}

test "curveBFromC releases allocations on failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkAllocationFailures, .{});
}

test "curveCFromB handles both roots and rejects singular and zero coefficients" {
    var p = try BigInt.initSet(std.testing.allocator, 41);
    defer p.deinit();
    var b = try BigInt.init(std.testing.allocator);
    defer b.deinit();
    for ([_]i32{ 14, 27, -14, 55 }) |value| {
        try b.set(value);
        var c = try curveCFromB(std.testing.allocator, &b, &p);
        defer c.deinit();
        try std.testing.expectEqual(@as(u32, 3), try c.toInt(u32));
    }
    for ([_]i32{ 0, 41 }) |value| {
        try b.set(value);
        try std.testing.expectError(error.InvalidCoefficient, curveCFromB(std.testing.allocator, &b, &p));
    }
    for ([_]i32{ 2, 39, -2 }) |value| {
        try b.set(value);
        try std.testing.expectError(error.SingularCurve, curveCFromB(std.testing.allocator, &b, &p));
    }
    try b.set(14);
    for ([_]i32{ -7, 0, 1, 2, 3, 4, 10 }) |value| {
        try p.set(value);
        try std.testing.expectError(error.InvalidPrime, curveCFromB(std.testing.allocator, &b, &p));
    }
}

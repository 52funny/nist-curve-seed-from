const std = @import("std");
const zeit = @import("zeit");
const Io = std.Io;

const nist_curve_seed_from = @import("nist_curve_seed_from");
pub const curveBFromC = @import("curve.zig").curveBFromC;
pub const curveCFromB = @import("curve.zig").curveCFromB;

const Drand48 = struct {
    x: u48 = 0,

    const a: u64 = 0x5deece66d;
    const c: u64 = 0xb;
    const m: u64 = 1 << 48;

    pub fn init(val: u64) Drand48 {
        const state = (val << 16) | 0x330e;
        return .{ .x = @truncate(state) };
    }

    pub fn next(self: *Drand48) u48 {
        // s[n+1] = (a * s[n] + b) % m
        const state = (a *% @as(u64, self.x) +% c) & (m - 1);
        self.x = @truncate(state);
        return self.x;
    }
};

pub fn main(init: std.process.Init) !void {
    const cases = [_]struct { curve_name: []const u8, time_str: []const u8, bits: usize, p: u521, b: u521 }{
        .{
            .curve_name = "P-256",
            .time_str = "1997-08-16 01:52:34",
            .bits = 256,
            .p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff,
            .b = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b,
        },
        .{
            .curve_name = "P-384",
            .time_str = "1997-08-30 14:50:08",
            .bits = 384,
            .p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff,
            .b = 0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef,
        },
    };

    const BigInt = std.math.big.int.Managed;

    for (cases) |case| {
        const timestamp = try parseUtcTimeStamp(case.time_str);
        var rng: Drand48 = .init(@bitCast(timestamp));
        var p = try BigInt.initSet(init.gpa, case.p);
        defer p.deinit();
        // hashSeed produces at most bits-1 bits, so c < p and direct equality
        // with the canonical target coefficient is sufficient.
        if (case.bits != p.bitCountAbs()) return error.InvalidBitLength;
        var expected_b = try BigInt.initSet(init.gpa, case.b);
        defer expected_b.deinit();
        var target_c = try curveCFromB(init.gpa, &expected_b, &p);
        defer target_c.deinit();

        var i: usize = 1;
        while (true) : (i += 1) {
            var seed_bytes: [20]u8 = undefined;
            for (&seed_bytes) |*byte| byte.* = @intCast(rng.next() >> 40);

            var c = try hashSeed(init.gpa, seed_bytes, case.bits);
            defer c.deinit();
            if (!c.eql(target_c)) continue;

            // Only derive a square root after a matching hash has been found.
            var b = curveBFromC(init.gpa, &c, &p) catch |err| switch (err) {
                error.NoSquareRoot, error.InvalidCoefficient, error.SingularCurve => continue,
                else => return err,
            };
            defer b.deinit();

            const actual = try b.toInt(u521);
            if (actual == case.b or (case.p % 4 == 1 and actual == case.p - case.b)) {
                std.debug.print("{s}: block={d}, timestamp={d}({s})\nseed=0x{x:0>40}, iters={}\nb=0x{x}\n\n", .{
                    case.curve_name,
                    i,
                    timestamp,
                    case.time_str,
                    std.mem.readInt(u160, &seed_bytes, .big),
                    i,
                    actual,
                });
                break;
            }
        }
    }
}

pub fn parseUtcTimeStamp(text: []const u8) !i64 {
    const inst = try zeit.instantFromText(.iso8601, text, &zeit.utc);
    return inst.unixTimestamp();
}

pub fn loopTimes(timestamp: u64, times: usize) ![20]u8 {
    if (times == 0) return error.InvalidIterationCount;
    var rng: Drand48 = .init(timestamp);
    var buff: [20]u8 = undefined;
    for (0..times) |_| {
        for (0..20) |i| {
            buff[i] = @intCast(rng.next() >> 40);
        }
    }
    return buff;
}

pub fn hashSeed(allocator: std.mem.Allocator, seed: [20]u8, l: usize) !std.math.big.int.Managed {
    if (l == 0) return error.InvalidBitLength;
    const v = (l - 1) / 160;
    const w = l - 160 * v - 1;

    const BigInt = std.math.big.int.Managed;
    const Sha1 = std.crypto.hash.Sha1;

    var digest: [20]u8 = undefined;
    Sha1.hash(&seed, &digest, .{});

    const mask = (@as(u160, 1) << @intCast(w)) - 1;

    var c = try BigInt.initSet(allocator, std.mem.readInt(u160, &digest, .big) & mask);
    errdefer c.deinit();

    var z = std.mem.readInt(u160, &seed, .big);
    var hi = try BigInt.init(allocator);
    defer hi.deinit();

    for (0..v) |_| {
        z +%= 1;

        var si: [20]u8 = undefined;
        std.mem.writeInt(u160, &si, z, .big);
        Sha1.hash(&si, &digest, .{});
        try hi.set(std.mem.readInt(u160, &digest, .big));

        // c = c || hi
        try c.shiftLeft(&c, 160);
        try c.bitOr(&c, &hi);
    }
    return c;
}

test "Drand48 should work" {
    const p256_date_time = "1997-08-16 01:52:34".*;
    const p256_timestamp = try parseUtcTimeStamp(&p256_date_time);

    try std.testing.expectEqual(@as(i64, 871696354), p256_timestamp);

    const buf = try loopTimes(@bitCast(p256_timestamp), 1842);
    const real_val = std.mem.readInt(u160, &buf, .big);
    try std.testing.expectEqual(0xC49D360886E704936A6678E1139D26B7819F7E90, real_val);
}

test "derive the published NIST curve coefficients from their seeds" {
    const cases = [_]struct { bits: usize, seed: u160, p: u521, b: u521 }{
        .{
            .bits = 192,
            .seed = 0x3045AE6FC8422F64ED579528D38120EAE12196D5,
            .p = 0xfffffffffffffffffffffffffffffffeffffffffffffffff,
            .b = 0x64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1,
        },
        .{
            .bits = 224,
            .seed = 0xBD71344799D5C7FCDC45B59FA3B9AB8F6A948BC5,
            .p = 0xffffffffffffffffffffffffffffffff000000000000000000000001,
            .b = 0xb4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4,
        },
        .{
            .bits = 256,
            .seed = 0xC49D360886E704936A6678E1139D26B7819F7E90,
            .p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff,
            .b = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b,
        },
        .{
            .bits = 384,
            .seed = 0xA335926AA319A27A1D00896A6773A4827ACDAC73,
            .p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffff0000000000000000ffffffff,
            .b = 0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef,
        },
        .{
            .bits = 521,
            .seed = 0xD09E8800291CB85396CC6717393284AAA0DA64BA,
            .p = std.math.maxInt(u521),
            .b = 0x0051953eb9618e1c9a1f929a21a0b68540eea2da725b99b315f3b8b489918ef109e156193951ec7e937b1652c0bd3bb1bf073573df883d2c34f1ef451fd46b503f00,
        },
    };
    for (cases) |case| {
        var seed: [20]u8 = undefined;
        std.mem.writeInt(u160, &seed, case.seed, .big);
        var c = try hashSeed(std.testing.allocator, seed, case.bits);
        defer c.deinit();
        var p = try std.math.big.int.Managed.initSet(std.testing.allocator, case.p);
        defer p.deinit();
        var b = try curveBFromC(std.testing.allocator, &c, &p);
        defer b.deinit();
        var published_b = try std.math.big.int.Managed.initSet(std.testing.allocator, case.b);
        defer published_b.deinit();
        var target_c = try curveCFromB(std.testing.allocator, &published_b, &p);
        defer target_c.deinit();
        try std.testing.expect(c.eql(target_c));
        const actual = try b.toInt(u521);
        // For P-224 either root is valid; p = 3 mod 4 uses the reference root.
        if (case.bits == 224) {
            try std.testing.expect(actual == case.b or actual == case.p - case.b);
        } else {
            try std.testing.expectEqual(case.b, actual);
        }
    }
}

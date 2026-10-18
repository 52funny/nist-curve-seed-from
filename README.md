# The NIST curves seed secret!


The original was posted on Anthropic blogs [https://www.anthropic.com/research/ecdsa-curve-origins](https://www.anthropic.com/research/ecdsa-curve-origins).


## Run Code!

For run the code, you should have Zig(0.16.0) enviroment.

Then, run it!

```bash
zig build run --Doptimize=ReleaseFast
```

The output is:

```text
P-256: block=1842, timestamp=871696354(1997-08-16 01:52:34)
seed=0xc49d360886e704936a6678e1139d26b7819f7e90, iters=1842
b=0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b

P-384: block=2688, timestamp=872952608(1997-08-30 14:50:08)
seed=0xa335926aa319a27a1d00896a6773a4827acdac73, iters=2688
b=0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef
```


## References

- [https://www.anthropic.com/research/ecdsa-curve-origins](https://www.anthropic.com/research/ecdsa-curve-origins)
- [https://words.filippo.io/seeds-solved](https://words.filippo.io/seeds-solved)
- [https://github.com/FiloSottile/mostly-harmless/blob/main/seeds/seeds.go](https://github.com/FiloSottile/mostly-harmless/blob/main/seeds/seeds.go)
- [https://std.neuromancer.sk/nist/P-256](https://std.neuromancer.sk/nist/P-256)
- [https://std.neuromancer.sk/nist/P-384](https://std.neuromancer.sk/nist/P-384)

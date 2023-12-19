import { BigInt } from './big-int';
import { consoleLog } from './env';

const _0 = BigInt.fromInt(0);
const _1 = BigInt.fromInt(1);
const _2 = BigInt.fromInt(2);
const _3 = BigInt.fromInt(3);
const _8 = BigInt.fromInt(8);
const POW_2_128 = _2.powInt(128);
const BETA = BigInt.fromString('55594575648329892869085402983802832744385952214688224221778511981742606582254');

export const P = BigInt.fromString('115792089237316195423570985008687907853269984665640564039457584007908834671663');
export const N = BigInt.fromString('115792089237316195423570985008687907852837564279074904382605163141518161494337');
export const GX = BigInt.fromString('55066263022277343669578718895168534326250603453777594175500187360389116729240');
export const GY = BigInt.fromString('32670510020758816978083085130507043184471273380659243275938904335757337482424');

export class Point {

    constructor(public readonly x: BigInt, public readonly y: BigInt) { }

    public static fromHex(hex: string): Point {
        const header = hex.substring(0, 2);

        if (hex.length !== 130 || header !== '04') {
            throw new Error('Invalid hex for point. Expected 65 bytes with the first byte beeing 0x04');
        }

        return new Point(
            BigInt.fromString(hex.substring(2, 66), 16),
            BigInt.fromString(hex.substring(66), 16)
        );
    }

    public multiply(k: BigInt): Point {
        const p = JPoint.fromAffine(this).multiplyWNAF(k);
    
        return p.toAffine();
        // return JPoint.fromAffine(this).multiplyWNAF(k).toAffine();
    }

    public add(p: Point): Point {
        return JPoint.fromAffine(this).add(JPoint.fromAffine(p)).toAffine();
    }

    public negate(): Point {
        return new Point(this.x, mod(this.y.clone().opposite()));
    }

    public equals(p: Point): boolean {
        return this.x === p.x && this.y === p.y;
    }

    public toHex(): string {
        return `04${ this.x.toString(16).padStart(64, '0') }${ this.y.toString(16).padStart(64, '0') }`;
    }

}

export const G: Point = new Point(GX, GY);

export const H = Point.fromHex('0450929b74c1a04954b78b4b6035e97a5e078a5a0f28ec96d547bfee9ace803ac031d3c6863973926e049e637cb1b5f40a36dac28af1766968c30c2313f3a38904');
export const B = Point.fromHex('043de7e317f561e8c9481b2128508c7effd2d524528b7da29e14d040d86e4b0159afd6259519fb77ba2b3bcb83a464cac183c85a2431539ad9a2ab41d7e06beeb2');

class JPoint {

    private static zero: JPoint = new JPoint(_0, _1, _0);

    private static base: JPoint = new JPoint(GX, GY, _1);

    private static instances: Map<BigInt, Map<BigInt, JPoint>> = new Map();

    private precomputedPoints: JPoint[] = [];

    constructor(public x: BigInt, public y: BigInt, public z: BigInt) { }

    public static fromAffine(p: Point): JPoint {
        if (p.x === _0 && p.y === _0) {
            return JPoint.zero;
        }

        if (p.x === GX && p.y === GY) {
            return JPoint.base;
        }

        const points: Map<BigInt, JPoint> = JPoint.instances.get(p.x) || new Map();

        let point = points.get(p.y);

        if (!point) {
            points.set(p.y, point = new JPoint(p.x, p.y, _1));

            JPoint.instances.set(p.x, points);
        }

        return point;
    }

    /**
     * Takes a bunch of Jacobian Points but executes only one
     * invert on all of them. invert is very slow operation,
     * so this improves performance massively.
     */
    static toAffineBatch(points: JPoint[]): Point[] {
        const toInv = invertBatch(points.map((p) => p.z));

        return points.map((p, i) => p.toAffine(toInv[i]));
    }

    static normalizeZ(points: JPoint[]): JPoint[] {
        return JPoint.toAffineBatch(points).map(JPoint.fromAffine);
    }

    /**
     * @see https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-1998-cmo-2
     */
    public add(p: JPoint): JPoint {
        const X1 = this.x;
        const Y1 = this.y;
        const Z1 = this.z;

        const X2 = p.x;
        const Y2 = p.y;
        const Z2 = p.z;

        if (X2 === _0 || Y2 === _0) {
            return this;
        }

        if (X1 === _0 || Y1 === _0) {
            return p;
        }

        // We're using same code in equals()
        const Z1Z1 = mod(Z1.clone().square());
        const Z2Z2 = mod(Z2.clone().square());
        const U1 = mod(X1.clone().mul(Z2Z2));
        const U2 = mod(X2.clone().mul(Z1Z1));
        const S1 = mod(mod(Y1.clone().mul(Z2)).mul(Z2Z2));
        const S2 = mod(mod(Y2.clone().mul(Z1)).mul(Z1Z1));
        const H = mod(U2.sub(U1));
        const r = mod(S2.sub(S1));

        // H = 0 meaning it's the same point.
        if (H === _0) {
            if (r === _0) {
                return this.double();
            }

            return JPoint.zero;
        }

        const HH = mod(H.square());
        const HHH = mod(H.clone().mul(HH));
        const V = mod(U1.mul(HH));
        const X3 = mod(r.clone().square().sub(HHH).sub(V.clone().mulInt(2)));
        const Y3 = mod(r.mul(V.sub(X3)).sub(S1.mul(HHH)));
        const Z3 = mod(Z1.clone().mul(Z2).mul(H));

        return new JPoint(X3, Y3, Z3);
    }

    /**
     * 
     * @see https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
     */
    public double(): JPoint {
        const X1 = this.x;
        const Y1 = this.y;
        const Z1 = this.z;

        const A = mod(X1.clone().square());
        const B = mod(Y1.clone().square());
        const C = mod(B.clone().square());
        const D = mod(mod(X1.clone().add(B).square()).sub(A).sub(C).mulInt(2));
        const E = mod(_3.clone().mul(A));
        const F = mod(E.clone().square());
        const X3 = mod(F.sub(D.clone().mulInt(2)));
        const Y3 = mod(E.mul(D.sub(X3)).sub(_8.clone().mul(C)));
        const Z3 = mod(Y1.clone().mul(Z1).mulInt(2));

        return new JPoint(X3, Y3, Z3);
    }

    public negate(): JPoint {
        return new JPoint(this.x, mod(this.y.clone().opposite()), this.z);
    }

    public multiply(k: BigInt): JPoint {
        if (k === _0) {
            return JPoint.zero;
        }

        if (k === _1) {
            return this;
        }

        // eslint-disable-next-line prefer-const
        const splitK = SplitScalar.split(k);

        let k1 = splitK.k1;
        let k2 = splitK.k2;

        let k1p = JPoint.zero;
        let k2p = JPoint.zero;

        // eslint-disable-next-line @typescript-eslint/no-this-alias
        let d: JPoint = this;

        while (k1 > _0 || k2 > _0) {
            if (k1.and(_1)) {
                k1p = k1p.add(d);
            }

            if (k2.and(_1)) {
                k2p = k2p.add(d);
            }

            d = d.double();
            k1.rightShift(1);
            k2.rightShift(1);
        }

        if (splitK.k1neg) {
            k1p = k1p.negate();
        }

        if (splitK.k2neg) {
            k2p = k2p.negate();
        }

        k2p = new JPoint(mod(k2p.x.mul(BETA)), k2p.y, k2p.z);

        return k1p.add(k2p);
    }

    public multiplyWNAF(k: BigInt, constantTime: boolean = false): JPoint {
        const splitK = SplitScalar.split(k);

        const points1 = this.wNAF(splitK.k1, constantTime);
        const points2 = this.wNAF(splitK.k2, constantTime);

        let k1p = points1[0];
        let k2p = points2[0];

        if (splitK.k1neg) {
            k1p = k1p.negate();
        }

        if (splitK.k2neg) {
            k2p = k2p.negate();
        }

        k2p = new JPoint(mod(k2p.x.mul(BETA)), k2p.y, k2p.z);

        // return JPoint.normalizeZ(constantTime ? [k1p.add(k2p), f1p.add(f2p)] : [k1p.add(k2p)])[0];

        constantTime && points1[1].add(points2[1]);

        return k1p.add(k2p);
    }

    public toAffine(invZ: BigInt = invert(this.z)): Point {
        const iz1 = invZ;
        const iz2 = mod(iz1.clone().square());
        const iz3 = mod(iz2.clone().mul(iz1));

        const ax = mod(this.x.clone().mul(iz2));
        const ay = mod(this.y.clone().mul(iz3));

        const zz = mod(this.z.clone().mul(iz1));

        if (zz !== _1) {
            throw new Error('invZ was invalid');
        }

        return new Point(ax, ay);
    }

    private wNAF(n: BigInt, constantTime: boolean = false): JPoint[] {
        // const W = this === JPoint.base ? 8 : 1;
        const W = 4;

        // Calculate precomputes on a first run, reuse them after
        const precomputes = this.getPrecomputedPoints(W);

        let p = JPoint.zero;
        let f = JPoint.zero;

        const windows = 1 + (128 / W); // W=8 17
        const windowSize = 2 ** (W - 1); // W=8 128
        const mask = BigInt.fromInt(2 ** W - 1); // Create mask with W ones: 0b11111111 for W=8
        const maxNumber = 2 ** W; // W=8 256

        for (let window = 0; window < windows; window++) {
            const offset = window * windowSize;

            // Extract W bits.
            let wbits = I32.parseInt(n.and(mask).toString());

            // Shift number by W bits.
            n.rightShift(W);

            // If the bits are bigger than max size, we'll split those.
            // +224 => 256 - 32
            if (wbits > windowSize) {
                wbits -= maxNumber;
                n.addInt(1);
            }

            // Check if we're onto Zero point.
            // Add random point inside current window to f.
            if (wbits === 0) {
                if (constantTime) {
                    // The most important part for const-time getPublicKey
                    let pr = precomputes[offset];

                    if (window % 2) {
                        pr = pr.negate();
                    }

                    f = f.add(pr);
                }
            } else {
                const cached = precomputes[offset + i32(Math.abs(wbits)) - 1];

                p = p.add(wbits < 0 ? cached.negate() : cached);
            }
        }

        return [p, f];
    }

    private getPrecomputedPoints(W: number): JPoint[] {
        if (this.precomputedPoints.length) {
            return this.precomputedPoints;
        }

        const windows = 128 / W + 1;
        const points: JPoint[] = this.precomputedPoints;

        // eslint-disable-next-line @typescript-eslint/no-this-alias
        let p: JPoint = this;
        let base: JPoint;

        for (let window = 0; window < windows; window++) {
            base = p;

            points.push(base);

            for (let i = 1; i < 2 ** (W - 1); i++) {
                base = base.add(p);

                points.push(base);
            }

            p = base.double();
        }

        // return this.precomputedPoints = JPoint.normalizeZ(points);
        return points;
    }

}

function mod(a: BigInt, b: BigInt = P): BigInt {
    const c = a.mod(b);

    return c >= _0 ? c : b.add(c);
}

function invert(aa: BigInt): BigInt {
    if (aa === _0) {
        throw new Error(`invert: expected positive integer, got ${ aa }`);
    }

    // Eucledian GCD https://brilliant.org/wiki/extended-euclidean-algorithm/
    let a = mod(aa.clone());

    let b = P;
    let x = _0;
    let y = _1;
    let u = _1;
    let v = _0;

    while (a !== _0) {
        const q = b.clone().div(a);
        const r = b.mod(a);
        const m = x.sub(u.mul(q));
        const n = y.sub(v.mul(q));

        b = a;
        a = r;
        x = u;
        y = v;
        u = m;
        v = n;
    }

    if (b !== _1) {
        throw new Error(`invert: does not exist, got ${ aa.toString() }`);
    }

    return mod(x);
}

const a1 = BigInt.fromString('3086d221a7d46bcde86c90e49284eb15', 16);
const b1 = BigInt.fromString('e4437ed6010e88286f547fa90abfe4c3', 16).opposite();
const a2 = BigInt.fromString('114ca50f7a8e2f3f657c1108d9d44cfd8', 16);
const b2 = a1;

class SplitScalar {

    constructor(
        public readonly k1: BigInt,
        public readonly k2: BigInt,
        public readonly k1neg: boolean,
        public readonly k2neg: boolean
    ) { }

    /** 
     * Split 256-bit K into 2 128-bit (k1, k2) for which k1 + k2 * lambda = K.
     * Used for endomorphism https://gist.github.com/paulmillr/eb670806793e84df628a7c434a873066
     */
    public static split(k: BigInt): SplitScalar {
        const c1 = divNearest(b2.clone().mul(k), N);
        const c2 = divNearest(b1.clone().mul(k).opposite(), N);

        let k1 = mod(k.clone().sub(c1.clone().mul(a1)).sub(c2.clone().mul(a2)), N);
        let k2 = mod(c1.mul(b1).opposite().sub(c2.mul(b2)), N);

        const k1neg = k1 > POW_2_128;
        const k2neg = k2 > POW_2_128;

        if (k1neg) {
            k1 = N.clone().sub(k1);
        }

        if (k2neg) {
            k2 = N.clone().sub(k2);
        }

        if (k1 > POW_2_128 || k2 > POW_2_128) {
            throw new Error('splitScalarEndo: Endomorphism failed, k=' + k.toString());
        }

        return new SplitScalar(k1, k2, k1neg, k2neg);
    }

}

function divNearest(a: BigInt, b: BigInt): BigInt {
    return a.add(b.clone().div(_2)).div(b);
}

/**
 * Takes a list of numbers, efficiently inverts all of them.
 * 
 * @param nums list of integers
 * 
 * @returns list of inverted integers
 */
function invertBatch(nums: BigInt[]): BigInt[] {
    const scratch: BigInt[] = new Array(nums.length);

    // Walk from first to last, multiply them by each other MOD p
    const lastMultiplied = nums.reduce(
        (acc, num, i) => {
            if (num === _0) {
                return acc;
            }

            scratch[i] = acc;

            return mod(acc.clone().mul(num));
        },
        _1.clone()
    );

    // Invert last element
    const inverted = invert(lastMultiplied);

    // Walk from last to first, multiply them by inverted each other MOD p
    nums.reduceRight(
        (acc, num, i) => {
            if (num === _0) {
                return acc;
            }

            scratch[i] = mod(acc.clone().mul(scratch[i]));

            return mod(acc.clone().mul(num));
        },
        inverted
    );

    return scratch;
}

// The entry file of your WebAssembly module.

import { GX, GY, Point } from './secp256k1';

export function multiply(p: string[], k: string): string[] {
    const q = new Point(BigInt.fromString(p[0], 16), BigInt.fromString(p[1], 16)).multiply(BigInt.fromString(k, 16));

    return [q.x.toString(16), q.y.toString(16)];
}

export const G = [GX.toString(16), GY.toString(16)];

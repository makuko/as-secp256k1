import benchmark from 'benchmark';
import BN from 'bn.js';
import elliptic from 'elliptic';

import { G, multiply } from '../build/release.js';

const suite = new benchmark.Suite('multiply', { minSamples: 1000});

const ec = new elliptic.ec('secp256k1');

// const f = random(N);
const f = BigInt(1234567890987654321n).toString(16);
const fBN = new BN(f, 16);

let oG = G.value;
let eG = ec.g;

suite
    .add('own', () => {
        oG = multiply(oG, f);
    })
    .add('elliptic', () => {
        eG = eG.mul(fBN);
    })
    .on('cycle', event => {
        console.log(String(event.target));
    })
    .on('complete', () => {
        console.log('Fastest is ' + suite.filter('fastest').map('name'));
    })
    .on('error', event => {
        console.error(event.target.error);
    })
    .run();

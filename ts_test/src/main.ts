const fs = require('fs');//내장 라이브러리 가져오기 include 비스무리한 거

const input = fs.readFileSync(0, 'utf8');//표준 입력 읽기
const lines = input.split(/\s+/);//입력 줄을 공백 등 기준으로 다 잘라서 배열로 저장

let cur = 0;

const N = parseInt(lines[0] || "0");
let _0 = 0;
let _1 = 0;
for (let i = 1; i <= N; i++) {
    const num = parseInt(lines[i] || "0");
    if (num === 0) {
        _0++;
    } else {
        _1++;
    }
}

console.log(Math.min(_0, _1));

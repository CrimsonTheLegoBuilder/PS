const fs = require('fs');//내장 라이브러리 가져오기 include 비스무리한 거

const input = fs.readFileSync(0, 'utf8');//표준 입력 읽기
const lines = input.split(/\s+/);//입력 줄을 공백 등 기준으로 다 잘라서 배열로 저장

let cur = 0;
while (cur < lines.length) {
    const nStr = lines[cur++];//하나 읽고 커서 내림
    if (!nStr) break;

    const N = parseInt(nStr);
    if (N === -1) break;

    let totalDistance = 0;
    let previousTime = 0;

    for (let i = 0; i < N; i++) {
        const s = parseInt(lines[cur++] || "0");//속도 읽기
        const t = parseInt(lines[cur++] || "0");//시간 읽기

        const duration = t - previousTime;//실제 주행 시간
        totalDistance += s * duration;
        previousTime = t;
    }

    process.stdout.write(totalDistance + " miles\n");
}

import re
import sys

def parse_point(s):
    # (x, y) 형태의 문자열을 파싱하여 좌표 반환
    match = re.search(r'\((-?\d+),\s*(-?\d+)\)', s)
    if match:
        return int(match.group(1)), int(match.group(2))
    return None

def check_sorted(label, cur, points, log_file):
    if not points or len(points) < 2:
        return
    
    # 캘리퍼스 방향(cur)에 수직인 벡터 (Normal Vector)
    # 이 축으로 투영했을 때 점들이 순서대로 정렬되어 있어야 함
    nx, ny = -cur[1], cur[0]
    
    projections = [p[0]*nx + p[1]*ny for p in points]
    
    # 오름차순 또는 내림차순인지 확인
    is_ascending = all(projections[i] <= projections[i+1] for i in range(len(projections)-1))
    is_descending = all(projections[i] >= projections[i+1] for i in range(len(projections)-1))
    
    # 정렬이 깨졌다면 로그 파일에 기록
    if not is_ascending and not is_descending:
        error_msg = (
            f"\n[ERROR] Sorting broken at {label}\n"
            f"Current Vector (cur): {cur}\n"
            f"Points: {points}\n"
            f"Projections: {projections}\n"
            f"{'-'*40}\n"
        )
        print(error_msg) # 콘솔 출력
        log_file.write(error_msg) # 파일 저장
        return False
    return True

def analyze_log(filename, output_filename):
    cur = None
    list_type = None
    points = []
    
    print(f"Analyzing {filename}...")
    
    try:
        with open(filename, 'r', encoding='utf-8') as f, open(output_filename, 'w', encoding='utf-8') as out_f:
            lines = f.readlines()
            
            for line in lines:
                line = line.strip()
                if not line: continue
                
                # 1. 현재 각도(cur) 갱신 시점
                if "DEBUG:: cur::" in line:
                    # 이전 리스트가 있었다면 정렬 검사 수행
                    if cur and points and list_type:
                        check_sorted(list_type, cur, points, out_f)
                    
                    cur = parse_point(line)
                    list_type = None
                    points = []
                    
                # 2. 리스트(bot.head, top.tail 등) 파싱
                elif "bot.head[" in line:
                    if list_type != "bot.head":
                        if cur and points: check_sorted(list_type, cur, points, out_f)
                        list_type = "bot.head"
                        points = []
                    p = parse_point(line)
                    if p: points.append(p)
                    
                elif "bot.tail[" in line:
                    if list_type != "bot.tail":
                        if cur and points: check_sorted(list_type, cur, points, out_f)
                        list_type = "bot.tail"
                        points = []
                    p = parse_point(line)
                    if p: points.append(p)
                    
                elif "top.head[" in line:
                    if list_type != "top.head":
                        if cur and points: check_sorted(list_type, cur, points, out_f)
                        list_type = "top.head"
                        points = []
                    p = parse_point(line)
                    if p: points.append(p)

                elif "top.tail[" in line:
                    if list_type != "top.tail":
                        if cur and points: check_sorted(list_type, cur, points, out_f)
                        list_type = "top.tail"
                        points = []
                    p = parse_point(line)
                    if p: points.append(p)

            # 마지막 블록 검사
            if cur and points and list_type:
                check_sorted(list_type, cur, points, out_f)
                
        print(f"Analysis Complete. Logs saved to {output_filename}")
        
    except FileNotFoundError:
        print(f"Error: Could not find file '{filename}'. Please check the file name.")

if __name__ == "__main__":
    # 분석할 파일명과 저장할 로그 파일명
    analyze_log('debug.txt', 'analysis_log.txt')
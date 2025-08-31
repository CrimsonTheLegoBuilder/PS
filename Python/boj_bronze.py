s=input().strip(); print(int(s,16) if s[:2]=='0x' else int(s,8) if s!='0' and s[0]=='0' else int(s))

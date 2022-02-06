## 양자내성 암호 알고리즘 구현

### ✔️ **프로젝트 개요**

- 타원곡선 기반 양자 내성 암호 구현
- 외부 라이브러리 사용을 최소화하고 클래스 구현

### 🗓️ 기간
2020년 가을학기 (수학과 캡스톤 프로젝트)
### 🏃🏻‍♂️ 진행 사항

- **소수 판별법**
    - Fermat, Miller-Rabin
- **인수 분해**
    - Pollard’s p-1 Factorization
    - Exponent Factoriztion
    - Quadratic Sieve
- **이산 로그**
    - BS-GS, CRT, Pohlig-Hellman
- **타원 곡선**
    - 타원 곡선 위의 점(Point), Field Element 간의 연산 정의
    - SHA-256 및 ECDSA 서명 알고리즘 구현
- **Polynomial Ring**
    - Polynomial Ring Element 및 원소들간의 연산 정의
- **NTRU**
    - 격자(lattice) 기반 양자 내성 암호(lattice problem)
    - Polynomial Ring을 통해 설계

### 💡 배운점 / 느낀점

- Class에 대한 깊은 이해
- 예외처리의 중요성
- 계산량 고려의 중요성
    - 큰 숫자간의 연산을 할 때 Numpy가 작동하지 않음
    - 시프트 연산을 활용해 문제 해결

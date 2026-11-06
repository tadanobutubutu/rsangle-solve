#!/usr/bin/env python3
"""rsangle (Daily AlpacaHack) の解法スクリプト。

使い方:
    python3 solve.py

配布ファイル output.txt に公開された n, e, c の値をスクリプト内に
含んでいます。外部への接続は不要です。
実行には Python 3.8 以降が必要です (pow(a, -1, m) で逆元を計算するため)。
因数分解には PARI/GP を使います。gp コマンドがあればそれを使い、
なければ Python バインディングの cypari2 を使います。どちらもない
環境では、SageMath の sage コマンドで因数分解します。いずれもない
場合は、pari-gp のインストール (例: apt install pari-gp) か、
pip install cypari2 か、SageMath のインストールが必要です。
"""

import re
import shutil
import subprocess
from math import gcd

# output.txt で公開された値。
# c = (c0, c1) は、剰余環 (Z/nZ)[x] / (x^2 + 1) の元 c0 + c1*i を表します。
# ここで i は x の剰余類で、i^2 = -1 を満たします。
n = 3824285698879494134357686342703215412032678268914632159283
e = 65537
c0 = 1656365636000604753217345256743048680199241083598544856864
c1 = 2614329834823602692469489358488649157330691937097719404460


# ---------------------------------------------------------------------------
# 剰余環 (Z/modZ)[i] での演算。元は (a, b) = a + b*i のペアで表します。
# ---------------------------------------------------------------------------

def ring_mul(left, right, modulus):
    a, b = left
    c, d = right
    return ((a * c - b * d) % modulus, (a * d + b * c) % modulus)


def ring_pow(value, exponent, modulus):
    result = (1, 0)
    while exponent:
        if exponent & 1:
            result = ring_mul(result, value, modulus)
        value = ring_mul(value, value, modulus)
        exponent >>= 1
    return result


def ring_inv(value, modulus):
    # (a + b*i)^(-1) = (a - b*i) / (a^2 + b^2)
    a, b = value
    inverse_norm = pow((a * a + b * b) % modulus, -1, modulus)
    return (a * inverse_norm % modulus, -b * inverse_norm % modulus)


# ---------------------------------------------------------------------------
# n の因数分解
# ---------------------------------------------------------------------------

def factor_with_sage(value):
    """SageMath の sage コマンドで value を因数分解し、素因数のリストを返します。"""
    if shutil.which("sage") is None:
        raise SystemExit(
            "因数分解には PARI/GP が必要です。pari-gp をインストールするか、"
            "pip install cypari2 を実行するか、SageMath をインストールしてください。"
        )
    result = subprocess.run(
        ["sage", "-c", f"print(factor({value}))"],
        text=True,
        capture_output=True,
        check=True,
    )
    # sage の factor の出力は "p * q" の形なので、数値を取り出します。
    return sorted(int(number) for number in re.findall(r"\d+", result.stdout))


def factor_n(value):
    """n を 2 つの素数に因数分解して (小さい方, 大きい方) を返します。"""
    if shutil.which("gp") is not None:
        # PARI/GP の factor で因数分解します。
        result = subprocess.run(
            ["gp", "-q", "-f"],
            input=f"factor({value})\n",
            text=True,
            capture_output=True,
            check=True,
        )
        rows = re.findall(r"\[\s*(\d+)\s+(\d+)\s*\]", result.stdout)
        primes = sorted(int(prime) for prime, exponent in rows if int(exponent) == 1)
    else:
        # gp がない環境では、PARI/GP の Python バインディングである
        # cypari2 を使います。どちらも同じ PARI の factor です。
        try:
            from cypari2 import Pari
        except ImportError:
            # cypari2 もない環境では、SageMath の sage コマンドで
            # 因数分解します。
            primes = factor_with_sage(value)
        else:
            pari = Pari()
            fac = pari.factor(pari(value))
            primes = sorted(
                int(fac[i, 0]) for i in range(fac.nrows()) if int(fac[i, 1]) == 1
            )

    assert len(primes) == 2 and primes[0] * primes[1] == value
    return primes[0], primes[1]


# ---------------------------------------------------------------------------
# 素数 mod での m の復元
# ---------------------------------------------------------------------------

def recover_m_mod_prime(prime, d):
    """c から z = c^d を復元し、m = (z - 1) / (i * (z + 1)) を取り出します。"""
    c_mod_prime = (c0 % prime, c1 % prime)

    # z は c の属するノルム 1 部分群の中で e 乗されており、
    # d はその部分群の位数に対する e の逆元なので、z = c^d で戻ります。
    z = ring_pow(c_mod_prime, d, prime)
    assert ring_pow(z, e, prime) == c_mod_prime

    # z = (1 + i*m) / (1 - i*m) を m について解いた式です。
    z_minus_one = ((z[0] - 1) % prime, z[1])
    z_plus_one = ((z[0] + 1) % prime, z[1])
    denominator = ring_mul((0, 1), z_plus_one, prime)  # i * (z + 1)
    m_pair = ring_mul(z_minus_one, ring_inv(denominator, prime), prime)

    # m は通常の整数なので、虚部は 0 になるはずです。
    assert m_pair[1] == 0
    return m_pair[0]


def main() -> None:
    p, q = factor_n(n)
    print(f"p = {p}")
    print(f"q = {q}")

    # z のノルムは 1 です。ノルム 1 の元全体は部分群をなし、その位数は
    # 素数ごとに決まります。素数を r として、x^2 + 1 が既約な場合
    # (r が mod 4 で 3) は r + 1、分解する場合 (r が mod 4 で 1) は r - 1 です。
    D = 1
    for prime in (p, q):
        D *= prime + 1 if prime % 4 == 3 else prime - 1
    print(f"norm-one subgroup order D = {D}")

    # e と D が互いに素なため、部分群の中では e 乗の逆演算ができます。
    assert gcd(e, D) == 1
    d = pow(e, -1, D)

    # 素数ごとに m を復元し、中国式剰余定理で合成します。
    m_p = recover_m_mod_prime(p, d)
    m_q = recover_m_mod_prime(q, d)
    m = m_p + p * ((m_q - m_p) * pow(p, -1, q) % q)
    assert 0 <= m < n

    flag = m.to_bytes((m.bit_length() + 7) // 8, "big").decode()
    print(f"flag = {flag}")


if __name__ == "__main__":
    main()

def circle_distance(a, b, N):
    d = abs(b - a)
    if d > N / 2:
        return N - d
    else:
        return d



def euler(fx, h, t0, tf):
    x = 0.
    for i in range((tf - t0) / h):
        x += x + h * fx(x)
    return x



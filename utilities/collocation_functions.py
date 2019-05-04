

# hermite interpolation
# ykpl := "y k plus 1"
def yc(yk, ykp1, fk, fkp1, hk):
    return .5 * (yk + ykp1) + (hk / 8.0) * (fk - fkp1)


def defect(yk, ykp1, fk, fkp1, fc, hk):
    return ykp1 - yk - (hk / 6.0) * (fk + 4.0 * fc + fkp1)


def state_vector_constructor(x, n, ns, gamma):
    y = []
    for i in range(gamma):
        y.append(x[n + i * (ns + 1) - 1])
    return y


# def center_state_vector_constructor(yk_full, ykp1_full, fk, fkp1, p, q, hk):
#     yc_array = []
#     for i in range(p):
#         yc_array.append(yc(yk_full[i], ykp1_full[i], fk, fkp1, hk))
#     for j in range(q):
#         yc_array.append((yk_full[p + j] + ykp1_full[p + j]) / 2)
#     return yc_array


def center_state_vector_constructor(yk_full, ykp1_full, fk, fkp1, p, q, hk, fxn_list):
    yc_array = []
    for i in range(p):
        if i % 2:
            idx = int(.5 * (i - 1))
            fxn = fxn_list[idx]
            yc_array.append(yc(yk_full[i], ykp1_full[i], fxn(yk_full), fxn(ykp1_full), hk))
        else:
            yc_array.append(0.0)
    for j in range(q):
        yc_array.append((yk_full[p + j] + ykp1_full[p + j]) / 2)
    return yc_array


def dfct_eval(x, fxn_list, n, ns, p, q, dfct_idx):
    """
    x: full NLP variable vector
    fxn: algebraic derivative expression
    n: segment number under evaluation
    ns: total number of segments
    p: number of states
    q: number of controls
    dfct_idx: defect index coresponding to specficic state
    """
    fxn = fxn_list[dfct_idx - 1]
    tf = x[-1]
    hk = tf / ns
    gamma = p + q
    yk_full = state_vector_constructor(x, n, ns, gamma)
    ykp1_full = state_vector_constructor(x, n + 1, ns, gamma)   
    # =================================================================    
    # asumes DAEs do not contain higher than second derivatives
    # further, index schema assumes all states have second derivative
    yk = yk_full[2 * dfct_idx - 2]
    ykp1 = ykp1_full[2 * dfct_idx - 2]
    # =================================================================   
    fk = fxn(yk_full)
    fkp1 = fxn(ykp1_full)
    ykc = center_state_vector_constructor(yk_full, ykp1_full, fk, fkp1, p, q, hk, fxn_list)
    fc = fxn(ykc)
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d


def constraint_list_gen(ns, p, q, fxn, constraint_fxn_list):
    # for use with python minimize format
    constraint_dict = {}
    c = 0
    for i in range(ns):
        for j in range(len(constraint_fxn_list)):
            c += 1
            args = (constraint_fxn_list, i + 1, ns, p, q, j + 1)
            # temp_dict = {'type': 'eq', 'fun': constraint_fxn_list[j], 'args': (i + 1,)}
            temp_dict = {'type': 'eq', 'fun': fxn, 'args': args}
            constraint_dict['con{}'.format(c)] = temp_dict

    return list(constraint_dict.values())


def bound_gen(bound_list, ns):
    bounds = []
    for i in range(len(bound_list)):
        for j in range(ns + 1):
            bounds.append(bound_list[i])
    return tuple(bounds)

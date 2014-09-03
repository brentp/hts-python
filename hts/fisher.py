from .htsffi import ffi, libhts

def fisher_exact_test(n11, n12, n21, n22):
    """
    >>> fisher_exact_test(20, 10, 10, 35) # doctest: +ELLIPSIS
    {'two_tail': 0.00023..., 'right_tail': 0.00014..., 'left_tail': 0.9999...}
    """
    left, right = ffi.new("double *"), ffi.new("double *")
    two = ffi.new("double *")
    libhts.kt_fisher_exact(n11, n12, n21, n22, left, right, two)
    return dict(left_tail=left[0], right_tail=right[0], two_tail=two[0])

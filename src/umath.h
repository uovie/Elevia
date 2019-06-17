/*
 *  Uovie Math
 */

namespace uovie {
    constexpr double pi = 3.141592653589793238;

    /*** factorial ***/
    int factorial(int val) {
        int res = 1;
        while (val > 1)
            res *= val--;
        return res;
    }

    /*** double factorial ***/
    int factorial2(int val) {
        int res = 1;
        for (int i = val; i >= 0; i = i - 2) {
            if (i == 0 || i == 1)
                break;
            else
                res *= i;
        }
        return res;
    }
}
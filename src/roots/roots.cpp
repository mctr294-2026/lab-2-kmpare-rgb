#include "roots.hpp"

// ------------------------------------------
// Object 1: Bisection Method
// ------------------------------------------
bool bisection(std::function<double(double)> f, double a, double b, double *root) {

    // Checks if there is a root
    if (f(a) * f(b) >= 0) {
        return false;
    }
    // Loop to find root
    for (int i = 0; i <1000000; i++) {

        // finds midpoint
        double c = (a + b) / 2.0;

        // updates the pointer to the best guess
        *root = c;

        double fc = f(c); // assigns variable to f(c)

        // if the absolute value of fc is close enough to zero it returns true
        if ((fc < 0 ? -fc : fc) < 1e-6) {
            return true;
        }

        // narrows down the interval
        if (f(a) * fc < 0) {
            b = c; // root in left half
        } else {
            a = c; // root in right half
        }
    }

    return false;
}

// ------------------------------------------
// Object 2: False Position Method
// ------------------------------------------
bool regula_falsi(std::function<double(double)> f, double a, double b, double *root) {

     // Checks if there is a root
    if (f(a) * f(b) >= 0) {
        return false;
    }
    // Loop to find root
    for (int i = 0; i <1000000; i++) {

        // assigns variables to make math cleaner
        double fa = f(a);
        double fb = f(b);

        // finds midpoint
        double c = a - (fa * (b - a)) / (fb - fa);

        // updates the pointer to the best guess
        *root = c;

        double fc = f(c); // assigns variable to f(c)

        // if the absolute value of fc is close enough to zero it returns true
        if ((fc < 0 ? -fc : fc) < 1e-6) {
            return true;
        }

        // narrows down the interval
        if (fa * fc < 0) {
            b = c; // root in left half
        } else {
            a = c; // root in right half
        }
    }
              
    return false;
}

// ------------------------------------------
// Object 3: Newton-Raphson Method
// ------------------------------------------
bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root) {
    
    double x = c;  // initial guess

    for (int i = 0; i < 1000000; i++) {

        double fx = f(x);
        double gx = g(x); // this is the derivitive
        
        // checks for derivitive of zero as this algorithm will fail
        if ((gx < 0 ? -gx : gx) < 1e-6) {
            return false;
        }

        // creates new guess
        double x_new = x - (fx / gx);

        *root = x_new;

        // checking if the change between the iterated guesses is within the tolerance
        double diff = x_new - x;
        if ((diff < 0 ? -diff : diff) < 1e-6) {
            return true;
        }

        x = x_new;
    }

    return false;
}

// ------------------------------------------
// Object 4: Secant Method
// ------------------------------------------
bool secant(std::function<double(double)> f, double a, double b, double c, double *root) {

    double x0 = a; // guess 1
    double x1 = b; // guess 2

    for (int i = 0; i < 1000000; i++) {

        double f0 = f(x0);
        double f1 = f(x1);

        // find the denominator
        double deno = f1 - f0;

        // double check to make sure deno is not zero
        if (( deno < 0 ? -deno : deno) < 1e-6) {
            return false;
        }

        // finds new guess
        double x_new = x1 - f1 * (x1 - x0) / deno;

        *root = x_new;

        // checking if the change between the iterated guesses is within the tolerance
        double diff = x_new - x1;
        if ((diff < 0 ? -diff : diff) < 1e-6) {
            return true;
        }

        x0 = x1;
        x1 = x_new;
    }
                
    return false;
}
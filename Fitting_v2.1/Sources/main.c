#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <C:/Users/lucas/CLionProjects/Calculator/tinyexpr.h>


// Funktion zum Berechnen und Ausgeben der Werte für eine gegebene Funktion y = f(x)
void berechne_funktion(const char *funktion, double x_start, double x_end, double schrittweite) {
    // Schleife über den Bereich von x und berechne f(x)
    for (double x = x_start; x <= x_end; x += schrittweite) {
        // Ersetze 'x' in der Funktion mit dem aktuellen Wert von x
        char ausdruck[256];
        snprintf(ausdruck, sizeof(ausdruck), funktion);

        // tinyExpr erwartet die Definition der Variablen explizit
        te_variable v[] = {{"x", &x}};  // Definiere x als Variable
        te_expr *expr = te_compile(ausdruck, v, 1, NULL); // kompiliere den Ausdruck mit der Variable x

        if (expr) {
            double ergebnis = te_eval(expr); // Berechne das Ergebnis
            printf("f(%.2f) = %.5f\n", x, ergebnis);
            te_free(expr); // Speicher freigeben
        } else {
            printf("Fehler beim Berechnen des Ausdrucks für x = %.2f\n", x);
        }
    }
}

int main() {
    char funktion[256];
    double x_start, x_end, schrittweite;

    // Benutzeraufforderung
    printf("Funktion f(x) (z.B. sin(x), x^2 + 3*x + 1): ");
    fgets(funktion, sizeof(funktion), stdin);

    funktion[strcspn(funktion, "\n")] = 0;

    // Eingabe des Bereichs und der Schrittweite
    printf("Bereich von x (Startwert und Endwert) (z.B. 0 10): ");
    scanf("%lf %lf", &x_start, &x_end);

    printf("Schrittweite (z.B. 0.1): ");
    scanf("%lf", &schrittweite);

    // Berechnung und Ausgabe der Funktionswerte
    berechne_funktion(funktion, x_start, x_end, schrittweite);

    return 0;
}

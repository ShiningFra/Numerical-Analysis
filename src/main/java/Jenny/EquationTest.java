package Jenny;

import org.junit.Assert;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class EquationTest {

    private static final double TOLERANCE = 1e-6; // Tolérance configurable

    public static void testResoudreEquation() {
        Imat solver = new EquationSolver();

        // Lire les cas de test depuis le fichier
        Integer z = 0;
        Integer y = 0;
        try (BufferedReader br = new BufferedReader(new FileReader("test_cases.txt"))) {
            String line;
            while ((line = br.readLine()) != null) {
                // Supposons que les lignes sont formatées comme "a b c"
                String[] parts = line.trim().split(" ");
                double a = getValueFromFeature(parts[0]);
                double b = getValueFromFeature(parts[1]);
                double c = getValueFromFeature(parts[2]);

                System.out.println("Test avec : a = " + a + ", b = " + b + ", c = " + c);
                
                // Résoudre l'équation
                try {
                    System.out.println("Début d'un cas de test");
                    y++;
                    double[] solutions = solver.resoudreEquation(a, b, c);

                    // Vérification des solutions
                    
                    for (double x : solutions) {
                        double result = a * x * x + b * x + c;
                        System.out.println("Résultat du test sur " + x + ": " + result);
                        Assert.assertEquals("L'équation ne donne pas 0 pour " + line + " avec x = " + x, 0, result, TOLERANCE);
                    }

                    // Si aucune solution n'est attendue, vérifier que le discriminant est négatif
                    if (solutions.length == 0) {
                        z++;
                        Assert.assertTrue("Le discriminant doit être négatif pour " + line, b * b - 4 * a * c < 0);
                    }
                    System.out.println("Fin d'un test");
                } catch (IllegalArgumentException e) {
                    System.out.println(e);
                }
            }
            System.out.println("Sucess rate " + (y - z) + " / " + y);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static double getValueFromFeature(String feature) {
        // Extraire la lettre et la convertir en valeur
        char letter = feature.charAt(0); // Par exemple, pour "a", on prend 'a'
        return Constants.getValue(letter);
    }

    public static void main(String[] args) {
        // Exécution de la résolution d'équation avec des valeurs d'exemple
        Imat solver = new EquationSolver();
        double a = 1; // Par exemple
        double b = -3;
        double c = 2;
        double[] result = solver.resoudreEquation(a, b, c);
        for (double r : result) {
            System.out.println("Solution : " + r);
        }
        System.out.println();
        testResoudreEquation();
    }
}
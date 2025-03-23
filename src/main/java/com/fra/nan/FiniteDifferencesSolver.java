/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.fra.nan;

import java.util.Arrays;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import javax.swing.JFrame;

public class FiniteDifferencesSolver {

    // Interface pour représenter une fonction réelle
    interface Function {
        double eval(double x);
    }
    
    public static void main(String[] args) {
        // Maillages à tester
        int[] N_values = {10, 20, 40, 80, 160, 320};

        // --- Cas 1 : u(x)=sin(pi*x) ---
        System.out.println("Cas 1: u(x)=sin(pi*x)");
        Function uExactSin = x -> Math.sin(Math.PI * x);
        Function fSin = x -> Math.PI * Math.PI * Math.sin(Math.PI * x) + Math.PI * Math.cos(Math.PI * x);
        double u0_sin = uExactSin.eval(0);
        double u1_sin = uExactSin.eval(1);
        runTests("sin(pi*x)", uExactSin, fSin, u0_sin, u1_sin, N_values);

        // --- Cas 2 : u(x)=x^3 ---
        System.out.println("Cas 2: u(x)=x^3");
        Function uExactCubic = x -> x * x * x;
        Function fCubic = x -> -6 * x + 3 * x * x;
        double u0_cubic = uExactCubic.eval(0);
        double u1_cubic = uExactCubic.eval(1);
        runTests("x^3", uExactCubic, fCubic, u0_cubic, u1_cubic, N_values);
    }

    /**
     * Pour un cas de test donné, résout l'équation sur différents maillages,
     * affiche l'erreur en norme Linfini et déduit l'ordre de convergence.
     * Pour le maillage le plus fin, affiche également un graphique comparant la solution exacte,
     * la solution approchée et l'erreur.
     */
    public static void runTests(String label, Function uExact, Function f, double u0, double u1, int[] N_values) {
        double prevError = -1;
        System.out.println("Test pour u(x)=" + label);
        double[] bestX = null, bestExact = null, bestApprox = null;
        int bestN = 0;
        for (int N : N_values) {
            double h = 1.0 / (N + 1);
            double[] x = new double[N + 2]; // incluant u0 et u1
            double[] u_ex = new double[N + 2];
            for (int i = 0; i < x.length; i++) {
                x[i] = i * h;
                u_ex[i] = uExact.eval(x[i]);
            }
            double[] u_app = solveFiniteDifferences(N, h, u0, u1, f);
            double error = computeMaxError(u_ex, u_app, N + 2);
            System.out.printf("  N = %d, h = %.5e, Erreur L_inf = %.5e", N, h, error);
            if (prevError > 0) {
                double order = Math.log(prevError / error) / Math.log(2);
                System.out.printf("  -> ordre ~ %.2f", order);
            }
            System.out.println();
            prevError = error;
            
            // Sauvegarde de la solution pour le maillage le plus fin
            bestX = x; 
            bestExact = u_ex; 
            bestApprox = u_app;
            bestN = N;
        }
        System.out.println("Affichage graphique pour N = " + bestN);
        plotGraph(bestX, bestExact, bestApprox, label);
    }

    /**
     * Résout l'équation -u'' + u' = f par différences finies pour un maillage de N points intérieurs.
     * La discrétisation utilisée est :
     *   (-1-h)*u_{i-1} + (2+h)*u_i - u_{i+1} = h^2 f(x_i), pour i = 1,...,N.
     * Les conditions aux limites u(0)=u0 et u(1)=u1 sont intégrées dans le second membre.
     */
    public static double[] solveFiniteDifferences(int N, double h, double u0, double u1, Function f) {
        double[] u = new double[N + 2]; // solution complète : u[0]=u0, u[N+1]=u1
        double[][] A = new double[N][N]; // matrice pour les points intérieurs
        double[] b = new double[N];       // second membre

        // Construction du système pour i = 1,...,N
        for (int i = 0; i < N; i++) {
            double x_i = (i + 1) * h;
            // Multiplication par h^2 de l'équation discrète
            b[i] = h * h * f.eval(x_i);
            // Coefficient de u_{i-1}
            if (i == 0) {
                // u0 connu
                b[i] -= (-1 - h) * u0;
            } else {
                A[i][i - 1] = -1 - h;
            }
            // Coefficient de u_i
            A[i][i] = 2 + h;
            // Coefficient de u_{i+1}
            if (i == N - 1) {
                // u1 connu
                b[i] -= (-1) * u1;
            } else {
                A[i][i + 1] = -1;
            }
        }

        // Résolution du système linéaire A * u_internal = b
        double[] u_internal = solveLinearSystem(A, b);
        // Assemblage de la solution complète
        u[0] = u0;
        System.arraycopy(u_internal, 0, u, 1, N);
        u[N + 1] = u1;
        return u;
    }

    /**
     * Résout le système linéaire A*x = b par élimination de Gauss (sans pivotage).
     */
    public static double[] solveLinearSystem(double[][] A, double[] b) {
        int N = b.length;
        double[] x = new double[N];
        // Élimination
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double factor = A[j][i] / A[i][i];
                for (int k = i; k < N; k++) {
                    A[j][k] -= factor * A[i][k];
                }
                b[j] -= factor * b[i];
            }
        }
        // Retour substitution
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }

    /**
     * Calcule l'erreur en norme L_infini (erreur maximale sur tous les points).
     */
    public static double computeMaxError(double[] u_exact, double[] u_approx, int len) {
        double maxError = 0.0;
        for (int i = 0; i < len; i++) {
            maxError = Math.max(maxError, Math.abs(u_exact[i] - u_approx[i]));
        }
        return maxError;
    }

    /**
     * Affiche un graphique avec trois courbes : solution exacte, solution approchée, et erreur.
     * Le graphique inclut une légende pour identifier chaque courbe.
     */
    public static void plotGraph(double[] x, double[] u_exact, double[] u_approx, String label) {
        XYSeries seriesExact = new XYSeries("Solution exacte");
        XYSeries seriesApprox = new XYSeries("Solution approchée");
        XYSeries seriesError = new XYSeries("Erreur");
        
        for (int i = 0; i < x.length; i++) {
            seriesExact.add(x[i], u_exact[i]);
            seriesApprox.add(x[i], u_approx[i]);
            seriesError.add(x[i], Math.abs(u_exact[i] - u_approx[i]));
        }
        
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(seriesExact);
        dataset.addSeries(seriesApprox);
        dataset.addSeries(seriesError);
        
        JFreeChart chart = ChartFactory.createXYLineChart(
                "Comparaison des solutions pour u(x) = " + label,
                "x", "u(x)",
                dataset, PlotOrientation.VERTICAL,
                true,   // Affichage de la légende
                true,   // Infobulles
                false   // URLs
        );
        
        // Forcer l'affichage de la légende (au cas où)
        if (chart.getLegend() != null) {
            chart.getLegend().setVisible(true);
        }
        
        // Création et affichage de la fenêtre graphique
        JFrame frame = new JFrame("Graphique pour u(x) = " + label);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.setSize(800, 600);  // Assure une taille suffisante pour voir la légende
        frame.setVisible(true);
    }
}

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.fra.nan;

/**
 *
 * @author Roddier
 */

import java.awt.BorderLayout;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Programme complet pour résoudre -u'' = f sur (0,1) par volumes finis,
 * estimer l'ordre de convergence en norme L∞ et afficher les graphiques.
 * On considère un maillage uniforme avec N mailles (donc, h = 1/N).
 */
public class FiniteVolumesSolver {

    /**
     * Résout le système discret obtenu avec la méthode des volumes finis.
     * Pour chaque cellule i, on considère le centre de contrôle à 
     * x_i = (i-0.5)*h pour i=1,...,N et on approxime les flux aux interfaces.
     * La formulation conduit à :
     * 
     *     u_{i-1} - 2 u_i + u_{i+1} = h^2 f(x_i)
     * 
     * pour i = 1,...,N, avec les conditions aux limites intégrées dans le second membre.
     *
     * @param N nombre de mailles (cellules)
     * @param u0 condition en x=0
     * @param u1 condition en x=1
     * @param f fonction f(x) (pour les valeurs au centre des cellules)
     * @return solution approchée aux centres de cellules, rangée de u1,...,u_N
     */
    public static double[] solveFV(int N, double u0, double u1, Functionn f) {
        double h = 1.0 / N;
        double[] x = new double[N]; // centres des cellules
        double[] b = new double[N];
        double[][] A = new double[N][N];

        // On considère le centre de la i-ème cellule: x_i = (i-0.5)*h, pour i = 1,...,N.
        // On assemble le système :
        //     u_{i-1} - 2 u_i + u_{i+1} = h^2 * f(x_i)
        // Avec u_0 = u(0) et u_{N+1} = u(1) connus.
        for (int i = 0; i < N; i++) {
            x[i] = (i + 0.5) * h; // centre de la cellule
            b[i] = - f.eval(x[i]) * h * h;
        }
        // Prise en compte des conditions limites
        b[0]    += u0;
        b[N - 1] += u1;

        // Construction de la matrice tridiagonale (schéma centré)
        for (int i = 0; i < N; i++) {
            A[i][i] = -2;
            if (i > 0) {
                A[i][i - 1] = 1;
            }
            if (i < N - 1) {
                A[i][i + 1] = 1;
            }
        }
        return gaussElimination(A, b);
    }

    /**
     * Méthode d'élimination de Gauss pour résoudre Ax = b.
     */
    public static double[] gaussElimination(double[][] A, double[] b) {
        int N = b.length;
        // Élimination
        for (int i = 0; i < N; i++) {
            for (int k = i + 1; k < N; k++) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < N; j++) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
        // Substitution arrière
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            x[i] = b[i];
            for (int j = i + 1; j < N; j++) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
        }
        return x;
    }

    /**
     * Calcule l'erreur en norme L∞ (erreur maximale) entre la solution exacte
     * évaluée aux centres de cellules et la solution approchée.
     *
     * @param uExact tableau de la solution exacte aux centres
     * @param uApprox tableau de la solution approchée aux centres
     * @return erreur L∞
     */
    public static double computeMaxError(double[] uExact, double[] uApprox) {
        double maxError = 0.0;
        for (int i = 0; i < uExact.length; i++) {
            maxError = Math.max(maxError, Math.abs(uExact[i] - uApprox[i]));
        }
        return maxError;
    }

    /**
     * Pour un cas de test donné, résout l'équation sur différents maillages,
     * affiche l'erreur en norme L∞ et estime l'ordre de convergence.
     * Pour le maillage le plus fin, affiche également un graphique comparant la solution exacte,
     * la solution approchée et l'erreur.
     *
     * @param label étiquette pour la fonction u(x)
     * @param uExact fonction u(x) exacte
     * @param uSecond fonction u''(x) (pour construire f(x) = -u''(x))
     * @param u0 condition à x=0
     * @param u1 condition à x=1
     * @param meshSizes tableau des nombres de mailles à tester
     */
    public static void runTests(String label, Functionn uExact, Functionn uSecond,
                                double u0, double u1, int[] meshSizes) {
        double prevError = -1;
        System.out.println("Test pour u(x) = " + uExact.getDescription());
        double[] bestX = null, bestExact = null, bestApprox = null;
        int bestN = 0;
        for (int N : meshSizes) {
            double h = 1.0 / N;
            // Pour les volumes finis, on évalue uExact aux centres de cellules.
            double[] x = new double[N];
            double[] uEx = new double[N];
            for (int i = 0; i < N; i++) {
                x[i] = (i + 0.5) * h;
                uEx[i] = uExact.eval(x[i]);
            }
            // f(x) = -u''(x)
            Functionn f = xVal -> -uSecond.eval(xVal);
            double[] uApp = solveFV(N, u0, u1, f);
            double error = computeMaxError(uEx, uApp);
            System.out.printf("N = %d, h = %.5e, erreur = %.5e", N, h, error);
            if (prevError > 0) {
                double order = Math.log(prevError / error) / Math.log(2);
                System.out.printf("  --> ordre ~ %.2f", order);
            }
            System.out.println();
            prevError = error;
            
            // On sauvegarde la solution pour le maillage le plus fin
            bestX = x; bestExact = uEx; bestApprox = uApp;
            bestN = N;
        }
        System.out.println("Affichage graphique pour N = " + bestN);
        plotGraph(bestX, bestExact, bestApprox, label);
    }

    /**
     * Point d'entrée principal.
     */
    public static void main(String[] args) {
        // Maillages à tester
        int[] meshSizes = {10, 20, 40, 80, 160, 320};

        // Cas 1 : u(x) = sin(pi*x)
        // u(0) = 0, u(1) = sin(pi)=0, et u''(x) = -pi^2 sin(pi*x)
        Functionn uExactSin = new Functionn() {
            public double eval(double x) {
                return Math.sin(Math.PI * x);
            }
            public String getDescription() {
                return "sin(pi*x)";
            }
        };
        Functionn uSecondSin = new Functionn() {
            public double eval(double x) {
                return -Math.PI * Math.PI * Math.sin(Math.PI * x);
            }
            public String getDescription() {
                return "u''(x) de sin(pi*x)";
            }
        };
        
        // Cas 2 : u(x) = x^3
        // u(0)=0, u(1)=1, et u''(x)=6x
        Functionn uExactCubic = new Functionn() {
            public double eval(double x) {
                return x * x * x;
            }
            public String getDescription() {
                return "x^3";
            }
        };
        Functionn uSecondCubic = new Functionn() {
            public double eval(double x) {
                return 6 * x;
            }
            public String getDescription() {
                return "u''(x) de x^3";
            }
        };

        System.out.println("u(x) = sin(pi*x)");
        runTests("sin(pi*x)", uExactSin, uSecondSin, uExactSin.eval(0), uExactSin.eval(1), meshSizes);
        System.out.println("u(x) = x^3");
        runTests("x^3", uExactCubic, uSecondCubic, uExactCubic.eval(0), uExactCubic.eval(1), meshSizes);
    }

    /**
     * Affiche un graphique avec trois courbes : solution exacte, solution approchée et erreur.
     * La légende identifie chacune des courbes.
     */
    public static void plotGraph(double[] x, double[] uExact, double[] uApprox, String label) {
        XYSeries seriesExact = new XYSeries("Solution exacte");
        XYSeries seriesApprox = new XYSeries("Solution approchée");
        XYSeries seriesError = new XYSeries("Erreur");
        for (int i = 0; i < x.length; i++) {
            seriesExact.add(x[i], uExact[i]);
            seriesApprox.add(x[i], uApprox[i]);
            seriesError.add(x[i], Math.abs(uExact[i] - uApprox[i]));
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(seriesExact);
        dataset.addSeries(seriesApprox);
        dataset.addSeries(seriesError);
        JFreeChart chart = ChartFactory.createXYLineChart(
                "Comparaison des solutions pour u(x) = " + label,
                "x", "u(x)",
                dataset,
                PlotOrientation.VERTICAL,
                true,  // légende
                true,  // infobulles
                false
        );
        showChart(chart, "Graphique pour u(x) = " + label);
    }

    /**
     * Affiche le graphique dans une fenêtre de taille réduite (600x400).
     */
    public static void showChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        ChartPanel chartPanel = new ChartPanel(chart);
        frame.getContentPane().add(chartPanel, BorderLayout.CENTER);
        frame.setSize(600, 400);
        frame.setVisible(true);
    }
}

/**
 * Interface pour définir une fonction réelle.
 */
interface Functionn {
    double eval(double x);
    default String getDescription() {
        return "";
    }
}

/**
 * Classe utilitaire pour afficher des graphiques à l'aide de JFreeChart.
 */
class Plotter {
    /**
     * Affiche deux séries (par exemple, la solution exacte et la solution approchée) dans un même graphique.
     * La légende identifie chacune des courbes.
     */
    public static void plotTwoSeries(double[] x, double[] series1, double[] series2,
                                     String title, String xLabel, String yLabel) {
        XYSeries s1 = new XYSeries("Exacte");
        XYSeries s2 = new XYSeries("Approchée");
        for (int i = 0; i < x.length; i++) {
            s1.add(x[i], series1[i]);
            s2.add(x[i], series2[i]);
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(s1);
        dataset.addSeries(s2);
        JFreeChart chart = ChartFactory.createXYLineChart(title, xLabel, yLabel,
                dataset, PlotOrientation.VERTICAL, true, true, false);
        showChart(chart, title);
    }

    /**
     * Affiche une seule série (par exemple, la courbe d'erreur) dans un graphique.
     * La légende permet d'identifier la série.
     */
    public static void plotSingleSeries(double[] x, double[] series, String title,
                                        String xLabel, String yLabel) {
        XYSeries s = new XYSeries(title);
        for (int i = 0; i < x.length; i++) {
            s.add(x[i], series[i]);
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(s);
        JFreeChart chart = ChartFactory.createXYLineChart(title, xLabel, yLabel,
                dataset, PlotOrientation.VERTICAL, true, true, false);
        showChart(chart, title);
    }

    /**
     * Affiche le graphique dans une fenêtre.
     * La taille de la fenêtre est fixée à 600x400.
     */
    private static void showChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        ChartPanel chartPanel = new ChartPanel(chart);
        frame.getContentPane().add(chartPanel, BorderLayout.CENTER);
        frame.setSize(600, 400);
        frame.setVisible(true);
    }
}


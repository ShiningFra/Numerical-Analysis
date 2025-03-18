/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.fra.nan;

/**
 *
 * @author Roddier
 */
import java.util.Arrays;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class FiniteDifferenceMethod {

    public static double fSin(double x) {
        return -Math.PI * Math.PI * Math.sin(Math.PI * x);
    }

    public static double[] solveBVP(int N, double u0, double u1, Function<Double, Double> f) {
        double h = 1.0 / (N + 1);
        double[] x = new double[N + 2];
        double[] u = new double[N];

        for (int i = 0; i < N + 2; i++) {
            x[i] = i * h;
        }

        double[][] A = new double[N][N];
        double[] b = new double[N];

        for (int i = 1; i <= N; i++) {
            A[i - 1][i - 1] = -2;
            if (i > 1) A[i - 1][i - 2] = 1;
            if (i < N) A[i - 1][i] = 1;
            b[i - 1] = f.apply(x[i]);
        }

        b[0] -= u0;
        b[N - 1] -= u1;

        double[] result = gaussElimination(A, b);
        return concatenate(u0, result, u1);
    }

    public static double[] gaussElimination(double[][] A, double[] b) {
        int n = b.length;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double factor = A[j][i] / A[i][i];
                for (int k = i; k < n; k++) {
                    A[j][k] -= factor * A[i][k];
                }
                b[j] -= factor * b[i];
            }
        }

        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            x[i] = b[i];
            for (int j = i + 1; j < n; j++) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
        }
        return x;
    }

    public static double[] concatenate(double u0, double[] u, double u1) {
        double[] result = new double[u.length + 2];
        result[0] = u0;
        System.arraycopy(u, 0, result, 1, u.length);
        result[result.length - 1] = u1;
        return result;
    }

    public static double computeError(double[] exact, double[] approx) {
        double maxError = 0.0;
        for (int i = 0; i < exact.length; i++) {
            maxError = Math.max(maxError, Math.abs(exact[i] - approx[i]));
        }
        return maxError;
    }

    public static void plotResults(double[] x, double[] exact, double[] approx) {
        XYSeries seriesExact = new XYSeries("Solution Exacte");
        XYSeries seriesApprox = new XYSeries("Solution Approchée");

        for (int i = 0; i < x.length; i++) {
            seriesExact.add(x[i], exact[i]);
            seriesApprox.add(x[i], approx[i]);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(seriesExact);
        dataset.addSeries(seriesApprox);

        JFreeChart chart = ChartFactory.createXYLineChart("Comparaison des Solutions", "x", "u(x)", dataset);
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame();
        frame.setContentPane(chartPanel);
        frame.pack();
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    public static void main(String[] args) {
        int[] meshSizes = {10, 20, 40, 80, 160, 320};
        double u0 = 0, u1 = 0; // Conditions aux limites
        double[] errors = new double[meshSizes.length];

        for (int i = 0; i < meshSizes.length; i++) {
            int N = meshSizes[i];
            double[] approx = solveBVP(N, u0, u1, FiniteDifferenceMethod::fSin);
            double[] exact = new double[N + 2];
            for (int j = 0; j < exact.length; j++) {
                exact[j] = Math.sin(Math.PI * j / (N + 1));
            }
            errors[i] = computeError(exact, approx);
            System.out.printf("Taille du maillage : %d, Erreur : %.6f%n", N, errors[i]);

            // Afficher le graphique pour le dernier maillage
            if (i == meshSizes.length - 1) {
                plotResults(exact, approx);
            }
        }

        // Calcul de l'ordre de convergence
        for (int i = 0; i < errors.length - 1; i++) {
            double order = Math.log(errors[i] / errors[i + 1]) / Math.log(meshSizes[i] / meshSizes[i + 1]);
            System.out.printf("Ordre de convergence entre %d et %d : %.6f%n", meshSizes[i], meshSizes[i + 1], order);
        }
    }
}

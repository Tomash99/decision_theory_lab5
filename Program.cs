using System;
using System.IO;
using System.Linq;
using System.Collections;
using System.Collections.Generic;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace matrixgame
{
    class Program
    {
        static double[,] filedata;

        static void Main(string[] args)
        {
            Console.WriteLine("Матрична гра");
            GetData();
            int index_min;
            int index_max;

            double max = Sidlo(out index_max);
            double min = Stovp(out index_min);
            if (max == min)
            {
                Console.WriteLine($"Ciдлова точка = {max} \nТип вирiшення: Чиста стратегiя!");
            }
            else
            {
                Console.WriteLine("Тип вирiшення: Змiшана стратегiя!");


                double[,] table = Zmishana();
                //for(int i = 0; i < 4; i++)
                //{
                //    for (int j = 0; j < 4; j++) 
                //    {
                //        Console.Write(table[i, j] + " ");
                //            }
                //    Console.WriteLine();
                //}

                double[] result = new double[3];
                double[,] table_result;
                Simplex S = new Simplex(table);
                table_result = S.Calculate(result);

                Console.WriteLine("Обрахована симплекс-таблиця:");
                for (int i = 0; i < table_result.GetLength(0); i++)
                {
                    for (int j = 0; j < table_result.GetLength(1); j++)
                        Console.Write(table_result[i, j] + " ");
                    Console.WriteLine();
                }

                Console.WriteLine();
                Console.WriteLine(" Оптимальний план :");
                Console.WriteLine("y[1] = " + result[0]);
                Console.WriteLine("y[2] = " + result[1]);
                Console.WriteLine("y[3] = " + result[2]);
                Console.WriteLine("Z(Y) = " + 1 * (result[0] + result[1] + result[2])+"\n");

                double[][] A = new double[][] { new double[] {1,4,5},
                                               new double[]{0,5,2},
                                               new double[]{0,3,11} };
                double[] vector = new double[3] { 0, 1, 1 }; 
                double[][] inv = MatrixInverse(A);
                Console.WriteLine();
                double[] VectorRESULT = new double[3];
                for (int i = 0; i < inv.Length; i++)
                {

                    VectorRESULT[i] =Math.Abs(vector[0] * inv[0][i] + vector[1] * inv[1][i] + vector[2] * inv[2][i]);
                }
                double Fx = VectorRESULT[0] + VectorRESULT[1] + VectorRESULT[2];
                Console.WriteLine($"x[1]={ VectorRESULT[0]} \nx[2]={ VectorRESULT[1]} \nx[3]={ VectorRESULT[2]}\n F(x)={Fx}");

                double p1 = Math.Abs( Fx * vector[0]);
                double p2 = Math.Abs(Fx * vector[1]);
                double p3 = Math.Abs( Fx * vector[2]);
                Console.WriteLine($"\nЦiна гри = {1 / Fx}");

                Console.WriteLine("\nОптимальна змiшана стратегiя для гравця 1");
                Console.WriteLine($"p1={p1} p2={p2} p3={p3}");

                double q1 = Math.Abs(Fx * result[0]);
                double q2 = Math.Abs( Fx * result[1]);
                double q3 = Math.Abs(Fx * result[2]);
                Console.WriteLine("\nОптимальна змiшана стратегiя для гравця 2");
                Console.WriteLine($"q1={q1} q2={q2} q3={q3}");
                Console.ReadLine();
            }

        }
        public static double[,] GetData()
        {
            string filePath = Path.GetFullPath("1.txt");
            int y = 0;
            int i = 0;
            using var fileReader = new StreamReader(filePath);
            string line;

            filedata = new double[5, 5];
            while ((line = fileReader.ReadLine()) != null)
            {
                string[] lines = line.Split(' ');


                for (i = 0; i < lines.Length; i++)
                {
                    filedata[y, i] = Convert.ToDouble(lines[i]);
                }
                y++;

            }
            return filedata;
        }

       
        static double Sidlo(out int index_max)
        {
            GetData();
            double x;
            double max;
            double[] mas;
            mas = new double[5];
            double[] mass;
            mass = new double[5];

            for (int y = 0; y < 5; y++)
            {
                for (int i = 0; i < 5; i++)
                {
                    mas[i] = filedata[y, i];
                }
                x = mas.Min();
                mass[y] = x;
            }
            max = mass.Max();
            index_max = Array.IndexOf(mass, max);
            return max;
        }

        static double Stovp(out int index_min)
        {
            GetData();
            double x;
            double min;
            double[] mas;
            mas = new double[5];
            double[] mass;
            mass = new double[5];

            for (int i = 0; i < 5; i++)
            {
                for (int y = 0; y < 5; y++)
                {
                    mas[y] = filedata[y, i];
                }
                x = mas.Max();
                mass[i] = x;
            }
            min = mass.Min();
            index_min = Array.IndexOf(mass, min);
            return min;
        }

        static double[,] Zmishana()
        {
            GetData();
            Kernel beforeChange, kernel = new Kernel(filedata);
            kernel.Print();
            do
            {
                beforeChange = kernel;
                for (int i = 0; i < filedata.GetLength(0); i++)
                {
                    for (int j = 0; j < filedata.GetLength(0); j++)
                    {
                        if (i == j)
                        {
                            continue;
                        }
                        if (kernel.StrategyMatrix[i] <= kernel.StrategyMatrix[j])
                        {
                            kernel.RemoveStrategy(i);
                        }
                    }
                }
                kernel = Kernel.Transpose(kernel);
                for (int i = 0; i < filedata.GetLength(1); i++)
                {
                    for (int j = 0; j < filedata.GetLength(1); j++)
                    {
                        if (i == j)
                        {
                            continue;
                        }

                        if (kernel.StrategyMatrix[i] >= kernel.StrategyMatrix[j])
                        {
                            kernel.RemoveStrategy(i);
                        }
                    }
                }
                kernel = Kernel.Transpose(kernel);
            } while (kernel != beforeChange);

            Console.WriteLine();

            Console.WriteLine();
            
            List<List<double>> numbers = new List<List<double>>() { };

            for (int i = 0; i < kernel.StrategyMatrix.Count; i++)
            {
                if (hasAppropriateNumber(kernel.StrategyMatrix[i].strategy))
                {
                    numbers.Add(new List<double>());
                }
                for (int j = 0; j < kernel.StrategyMatrix[i].strategy.Count; j++)
                {
                    if (kernel.StrategyMatrix[i].strategy[j] < 1000)
                    {
                        numbers[numbers.Count - 1].Add(kernel.StrategyMatrix[i].strategy[j]);
                        Console.Write(numbers[numbers.Count - 1][numbers[numbers.Count - 1].Count - 1]+" ");

                    }
                                    }
                Console.WriteLine();
            }

            double[,] matr = new double[numbers.Count + 1, numbers[0].Count + 1];
            for (int i = 0; i < (numbers.Count) + 1; i++)
            {
                if (i == 3) { matr[i, 0] = 0; }
                else { matr[i, 0] = 1; }

                for (int j = 1; j < (numbers.Count) + 1; j++)
                {
                    if (i == 3) { matr[i, j] = -1; }
                    else { matr[i, j] = numbers[i][j - 1]; }
                    
                }
            }
            return matr;
        }
        public static bool hasAppropriateNumber(List<double> arr)
        {
            for (int i = 0; i < arr.Count; i++)
            {
                if (arr[i] < 1000)
                {
                    return true;
                }
            }
            return false;
        }
        static double[][] CopyMatrix(double[,] matrix)
        {
            double[][] CopyA=new double[][] { };
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    CopyA[i][j] = matrix[i,j];
                }
            }
            return CopyA;
        }
        static double[][] MatrixInverse(double[][] matrix)
        {
            int n = matrix.Length;
            double[][] result = MatrixDuplicate(matrix);

            int[] perm;
            int toggle;
            double[][] lum = MatrixDecompose(matrix, out perm,
              out toggle);
            if (lum == null)
                throw new Exception("Unable to compute inverse");

            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;
                }

                double[] x = HelperSolve(lum, b);

                for (int j = 0; j < n; ++j)
                    result[j][i] = x[j];
            }
            return result;
        }

        static double[][] MatrixDuplicate(double[][] matrix)
        {
            double[][] result = matrix;
            for (int i = 0; i < matrix.Length; ++i) 
                for (int j = 0; j < matrix[i].Length; ++j)
                    result[i][j] = matrix[i][j];
            return result;
        }

        static double[] HelperSolve(double[][] luMatrix, double[] b)
        {
            int n = luMatrix.Length;
            double[] x = new double[n];
            b.CopyTo(x, 0);

            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum;
            }

            x[n - 1] /= luMatrix[n - 1][n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum / luMatrix[i][i];
            }

            return x;
        }

        static double[][] MatrixDecompose(double[][] matrix, out int[] perm, out int toggle)
        {
            int rows = matrix.Length;
            int cols = matrix[0].Length; 
            if (rows != cols)
                throw new Exception("Attempt to decompose a non-square m");

            int n = rows; 

            double[][] result = MatrixDuplicate(matrix);

            perm = new int[n]; 
            for (int i = 0; i < n; ++i) { perm[i] = i; }

            toggle = 1; 
            for (int j = 0; j < n - 1; ++j) 
            {
                double colMax = Math.Abs(result[j][j]); 
                int pRow = j;
         
                for (int i = j + 1; i < n; ++i)
                {
                    if (Math.Abs(result[i][j]) > colMax)
                    {
                        colMax = Math.Abs(result[i][j]);
                        pRow = i;
                    }
                }

                if (pRow != j)
                {
                    double[] rowPtr = result[pRow];
                    result[pRow] = result[j];
                    result[j] = rowPtr;

                    int tmp = perm[pRow]; 
                    perm[pRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle; 
                }

                if (result[j][j] == 0.0)
                {
                    int goodRow = -1;
                    for (int row = j + 1; row < n; ++row)
                    {
                        if (result[row][j] != 0.0)
                            goodRow = row;
                    }

                    if (goodRow == -1)
                        throw new Exception("Cannot use Doolittle's method");

                    double[] rowPtr = result[goodRow];
                    result[goodRow] = result[j];
                    result[j] = rowPtr;

                    int tmp = perm[goodRow]; 
                    perm[goodRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle; 
                }
               
                for (int i = j + 1; i < n; ++i)
                {
                    result[i][j] /= result[j][j];
                    for (int k = j + 1; k < n; ++k)
                    {
                        result[i][k] -= result[i][j] * result[j][k];
                    }
                }
            }
            return result;
        }
    }
}



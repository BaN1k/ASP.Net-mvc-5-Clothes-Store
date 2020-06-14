using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Numerics;

namespace AmcLab0
{

    interface IObgect
    {
        void GetPoints();
        void GetTriangles();
        void WriteFile(string path);
        void GetSplit(double nxa, double nya, double nza, double nxb, double nyb, double nzb);
        void GetSquare();
        void WriteFileSplit(string path1, string path2);

    }

    public class Cube : IObgect
    {
        //1drawing
        public const int n = 8, m = 3;
        public double x, y, z, xp, yp, zp, l;
        public double[,] points = new double[n, m];
        public double[,] triangles = new double[36, m];
        //2spliting
        public int todo = 1, q = 0, p = 0, itr1 = 0, itr2 = 0;
        public double[,] pointsplit1 = new double[36, m];
        public double[,] pointsplit;
        public double[,] pointsall;
        public double[,] trisp1;
        public double[,] trisp2;
        public double[,] sqr;
        public double nxa, nya, nza, nxb, nyb, nzb, A, B, C, D, T, x1, y1, z1, L = 0, M = 0, N = 0, mx, my, mz, l1, l2 = 0, l3 = 0;
        private int leng;

        public double[,] normals = {
            {-1,0,0 },
            {1,0,0 },
            {0,-1,0 },
            {0,1,0 },
            {0,0,-1 },
            {0,0,1 },
        };
        public Cube(double l, double x, double y, double z)
        {
            this.l = l;
            this.x = x;
            this.y = y;
            this.z = z;
            xp = l + x;
            yp = l + y;
            zp = l + z;
        }
        public void GetPoints()
        {
            int p = 0;
            for (double i = x; i <= xp; i += l)
            {
                for (double j = y; j <= yp; j += l)
                {
                    for (double k = z; k <= zp; k += l)
                    {
                        points[p, 0] = i;
                        points[p, 1] = j;
                        points[p, 2] = k;
                        p++;
                    }
                }
            }
        }
        public void GetTriangles()
        {
            double[,] dots = { { x, xp }, { y, yp }, { z, zp } };
            int p1 = 0, p2 = 0, p = 1;
            for (int j = 0; j < 3; j++)
            {
                p1 = p2;
                p2 += 6;
                for (int i = 0; i < n; i++)
                {
                    if (points[i, j] == dots[j, 0])
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            triangles[p1, k] = points[i, k];
                        }
                        p1++;
                    }
                    if (points[i, j] == dots[j, 1])
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            triangles[p2, k] = points[i, k];
                        }
                        p2++;
                    }
                }
                p2 += 2;
            }
            for (int i = 0; i < 36; i++)
            {
                if (p % 4 == 0)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        triangles[i + 2, k] = triangles[i, k];
                        triangles[i + 1, k] = triangles[i - 1, k];
                        triangles[i, k] = triangles[i - 2, k];
                    }
                }
                p++;
                if ((i + 1) % 6 == 0) p = 1;
            }
        }
        public void WriteFile(string path)
        {
            int p = 1, pn = 0;
            using (StreamWriter sw = new StreamWriter(path))
            {
                sw.WriteLine("solid Cube");
                for (int i = 0; i < 36; i++)
                {
                    if (p == 1)
                    {
                        sw.Write("facet normal");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write(Convert.ToString("\t" + normals[pn, j].ToString(CultureInfo.InvariantCulture)));
                        }
                        sw.WriteLine();
                        sw.WriteLine("outer loop");
                    }
                    sw.Write("vertex");
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(Convert.ToString("\t" + triangles[i, j].ToString(CultureInfo.InvariantCulture)));
                    }
                    sw.Write("\n");
                    if (p == 3)
                    {
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }
                    p++;
                    if ((i + 1) % 3 == 0) p = 1;
                    if ((i + 1) % 6 == 0) pn++;
                }
                sw.WriteLine("endsolid Cube");
                sw.Close();
            }
            Console.WriteLine("File succes created");
        }
        public void GetSplit(double _nxa, double _nya, double _nza, double _nxb, double _nyb, double _nzb)
        {
            nxa = _nxa;
            nya = _nya;
            nza = _nza;
            nxb = _nxb;
            nyb = _nyb;
            nzb = _nzb;
            int it1 = 0, it2 = 0, it3 = 0, it4 = 0, chk = 0;
            double vx, vy, vz, nxc, nyc, nzc, length, lb1, lc1, lb2, lc2, b1, b2, b3, c1, c2, c3;
            double[,] temp1 = new double[36, 3];
            double[,] temp2 = new double[36, 3];
            double[,] temp3 = new double[36, 3];
            double[,] temp4 = new double[36, 3];
            vx = nxb - nxa;
            vy = nyb - nya;
            vz = nzb - nza;
            length = Math.Sqrt(Math.Pow((vx), 2) + Math.Pow((vy), 2) + Math.Pow((vz), 2));
            vx = vx / length;
            vy = vy / length;
            vz = vz / length;
            nxc = nxa - vx * length;
            nyc = nya - vy * length;
            nzc = nza - vz * length;
            pointsall = this.triangles;
            A = nxb - nxa;
            B = nyb - nya;
            C = nzb - nza;
            D = A * (-nxa) + B * (-nya) + C * (-nza);
            for (int i = 0; i < 36; i += 3)
            {
                chk = 0;
                b1 = Math.Sqrt(Math.Pow((nxb - pointsall[i, 0]), 2) + Math.Pow((nyb - pointsall[i, 1]), 2) + Math.Pow((nzb - pointsall[i, 2]), 2));
                c1 = Math.Sqrt(Math.Pow((nxc - pointsall[i, 0]), 2) + Math.Pow((nyc - pointsall[i, 1]), 2) + Math.Pow((nzc - pointsall[i, 2]), 2));
                b2 = Math.Sqrt(Math.Pow((nxb - pointsall[i + 1, 0]), 2) + Math.Pow((nyb - pointsall[i + 1, 1]), 2) + Math.Pow((nzb - pointsall[i + 1, 2]), 2));
                c2 = Math.Sqrt(Math.Pow((nxc - pointsall[i + 1, 0]), 2) + Math.Pow((nyc - pointsall[i + 1, 1]), 2) + Math.Pow((nzc - pointsall[i + 1, 2]), 2));
                b3 = Math.Sqrt(Math.Pow((nxb - pointsall[i + 2, 0]), 2) + Math.Pow((nyb - pointsall[i + 2, 1]), 2) + Math.Pow((nzb - pointsall[i + 2, 2]), 2));
                c3 = Math.Sqrt(Math.Pow((nxc - pointsall[i + 2, 0]), 2) + Math.Pow((nyc - pointsall[i + 2, 1]), 2) + Math.Pow((nzc - pointsall[i + 2, 2]), 2));
                if (b1 < c1 && b2 < c2 && b3 < c3)
                {
                    temp3[it3, 0] = pointsall[i, 0];
                    temp3[it3, 1] = pointsall[i, 1];
                    temp3[it3, 2] = pointsall[i, 2];
                    it3++;
                    temp3[it3, 0] = pointsall[i + 1, 0];
                    temp3[it3, 1] = pointsall[i + 1, 1];
                    temp3[it3, 2] = pointsall[i + 1, 2];
                    it3++;
                    temp3[it3, 0] = pointsall[i + 2, 0];
                    temp3[it3, 1] = pointsall[i + 2, 1];
                    temp3[it3, 2] = pointsall[i + 2, 2];
                    it3++;
                }
                else if (b1 > c1 && b2 > c2 && b3 > c3)
                {
                    temp4[it4, 0] = pointsall[i, 0];
                    temp4[it4, 1] = pointsall[i, 1];
                    temp4[it4, 2] = pointsall[i, 2];
                    it4++;
                    temp4[it4, 0] = pointsall[i + 1, 0];
                    temp4[it4, 1] = pointsall[i + 1, 1];
                    temp4[it4, 2] = pointsall[i + 1, 2];
                    it4++;
                    temp4[it4, 0] = pointsall[i + 2, 0];
                    temp4[it4, 1] = pointsall[i + 2, 1];
                    temp4[it4, 2] = pointsall[i + 2, 2];
                    it4++;
                }

                for (int j = 0; j < 2; j++)
                {
                    for (int k = 1; k < 3; k++)
                    {
                        if (j != k)
                        {
                            L = pointsall[i + j, 0] - pointsall[i + k, 0];
                            M = pointsall[i + j, 1] - pointsall[i + k, 1];
                            N = pointsall[i + j, 2] - pointsall[i + k, 2];
                            l1 = Math.Sqrt(Math.Pow(L, 2) + Math.Pow(M, 2) + Math.Pow(N, 2));
                            x1 = pointsall[i + j, 0];
                            y1 = pointsall[i + j, 1];
                            z1 = pointsall[i + j, 2];
                            if ((A * L + B * M + C * N) != 0)
                            {
                                T = -((A * x1 + B * y1 + C * z1 + D) / (A * L + B * M + C * N));
                                mx = L * T + x1;
                                my = M * T + y1;
                                mz = N * T + z1;
                                l2 = Math.Sqrt(Math.Pow((pointsall[i + j, 0] - mx), 2) + Math.Pow((pointsall[i + j, 1] - my), 2) + Math.Pow((pointsall[i + j, 2] - mz), 2));
                                l3 = Math.Sqrt(Math.Pow((pointsall[i + k, 0] - mx), 2) + Math.Pow((pointsall[i + k, 1] - my), 2) + Math.Pow((pointsall[i + k, 2] - mz), 2));
                                todo = 1;
                                for (int e = 0; e < q + 1; e++) { if (pointsplit1[e, 0] == mx && pointsplit1[e, 1] == my && pointsplit1[e, 2] == mz) todo = 0; }
                                bool bi = false;
                                double dep = l1 - (l2 + l3);
                                if (dep > 0 && dep < 0.001) { bi = true; }
                                else if (dep < 0 && dep > -0.001) { bi = true; }
                                else if (dep == 0) { bi = true; }
                                if (bi && todo == 1)
                                {
                                    pointsplit1[q, 0] = mx;
                                    pointsplit1[q, 1] = my;
                                    pointsplit1[q, 2] = mz;
                                    ++q;
                                }
                                if (bi)
                                {
                                    lb1 = Math.Sqrt(Math.Pow((nxb - pointsall[i + j, 0]), 2) + Math.Pow((nyb - pointsall[i + j, 1]), 2) + Math.Pow((nzb - pointsall[i + j, 2]), 2));
                                    lc1 = Math.Sqrt(Math.Pow((nxc - pointsall[i + j, 0]), 2) + Math.Pow((nyc - pointsall[i + j, 1]), 2) + Math.Pow((nzc - pointsall[i + j, 2]), 2));
                                    lb2 = Math.Sqrt(Math.Pow((nxb - pointsall[i + k, 0]), 2) + Math.Pow((nyb - pointsall[i + k, 1]), 2) + Math.Pow((nzb - pointsall[i + k, 2]), 2));
                                    lc2 = Math.Sqrt(Math.Pow((nxc - pointsall[i + k, 0]), 2) + Math.Pow((nyc - pointsall[i + k, 1]), 2) + Math.Pow((nzc - pointsall[i + k, 2]), 2));
                                    if ((lb1 < lc1 && lb2 > lc2) || (lb1 < lc1 && lb2 == lc2))
                                    {
                                        if (chk == 0)
                                        {
                                            //1
                                            temp1[it1, 0] = pointsall[i + j, 0];
                                            temp1[it1, 1] = pointsall[i + j, 1];
                                            temp1[it1, 2] = pointsall[i + j, 2];
                                            ++it1;
                                            temp1[it1, 0] = mx;
                                            temp1[it1, 1] = my;
                                            temp1[it1, 2] = mz;
                                            ++it1;
                                            //2
                                            temp2[it2, 0] = pointsall[i + k, 0];
                                            temp2[it2, 1] = pointsall[i + k, 1];
                                            temp2[it2, 2] = pointsall[i + k, 2];
                                            ++it2;
                                            temp2[it2, 0] = mx;
                                            temp2[it2, 1] = my;
                                            temp2[it2, 2] = mz;
                                            ++it2;
                                            chk++;

                                        }
                                        else if (chk == 1)
                                        {
                                            if (temp1[it1 - 2, 0] == pointsall[i + j, 0] && temp1[it1 - 2, 1] == pointsall[i + j, 1] && temp1[it1 - 2, 2] == pointsall[i + j, 2])
                                            {
                                                //1
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                //2
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                                temp2[it2, 0] = temp2[it2 - 3, 0];
                                                temp2[it2, 1] = temp2[it2 - 3, 1];
                                                temp2[it2, 2] = temp2[it2 - 3, 2];
                                                ++it2;
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                                temp2[it2, 0] = pointsall[i + k, 0];
                                                temp2[it2, 1] = pointsall[i + k, 1];
                                                temp2[it2, 2] = pointsall[i + k, 2];
                                                ++it2;
                                            }
                                            else
                                            {
                                                //1
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                temp1[it1, 0] = temp1[it1 - 3, 0];
                                                temp1[it1, 1] = temp1[it1 - 3, 1];
                                                temp1[it1, 2] = temp1[it1 - 3, 2];
                                                ++it1;
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                temp1[it1, 0] = pointsall[i + j, 0];
                                                temp1[it1, 1] = pointsall[i + j, 1];
                                                temp1[it1, 2] = pointsall[i + j, 2];
                                                ++it1;
                                                //2
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                            }
                                        }
                                    }
                                    else if ((lb1 > lc1 && lb2 < lc2) || (lb1 == lc1 && lb2 < lc2))
                                    {
                                        if (chk == 0)
                                        {
                                            //1
                                            temp1[it1, 0] = pointsall[i + k, 0];
                                            temp1[it1, 1] = pointsall[i + k, 1];
                                            temp1[it1, 2] = pointsall[i + k, 2];
                                            ++it1;
                                            temp1[it1, 0] = mx;
                                            temp1[it1, 1] = my;
                                            temp1[it1, 2] = mz;
                                            ++it1;
                                            //2
                                            temp2[it2, 0] = pointsall[i + j, 0];
                                            temp2[it2, 1] = pointsall[i + j, 1];
                                            temp2[it2, 2] = pointsall[i + j, 2];
                                            ++it2;
                                            temp2[it2, 0] = mx;
                                            temp2[it2, 1] = my;
                                            temp2[it2, 2] = mz;
                                            ++it2;
                                            chk++;

                                        }
                                        else if (chk == 1)
                                        {
                                            if (temp1[it1 - 2, 0] == pointsall[i + k, 0] && temp1[it1 - 2, 1] == pointsall[i + k, 1] && temp1[it1 - 2, 2] == pointsall[i + k, 2])
                                            {
                                                //1
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                //2
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                                temp2[it2, 0] = temp2[it2 - 3, 0];
                                                temp2[it2, 1] = temp2[it2 - 3, 1];
                                                temp2[it2, 2] = temp2[it2 - 3, 2];
                                                ++it2;
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                                temp2[it2, 0] = pointsall[i + j, 0];
                                                temp2[it2, 1] = pointsall[i + j, 1];
                                                temp2[it2, 2] = pointsall[i + j, 2];
                                                ++it2;
                                            }
                                            else
                                            {
                                                //1
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                temp1[it1, 0] = temp1[it1 - 3, 0];
                                                temp1[it1, 1] = temp1[it1 - 3, 1];
                                                temp1[it1, 2] = temp1[it1 - 3, 2];
                                                ++it1;
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                temp1[it1, 0] = pointsall[i + k, 0];
                                                temp1[it1, 1] = pointsall[i + k, 1];
                                                temp1[it1, 2] = pointsall[i + k, 2];
                                                ++it1;
                                                //2
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            itr1 = it1 + it3;
            itr2 = it2 + it4;
            trisp1 = new double[itr1, m];
            trisp2 = new double[itr2, m];

            for (int i = 0; i < it1; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    trisp1[i, j] = temp1[i, j];
                }
            }
            for (int i = 0; i < it3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    trisp1[i + it1, j] = temp3[i, j];
                }
            }
            for (int i = 0; i < it2; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    trisp2[i, j] = temp2[i, j];
                }
            }
            for (int i = 0; i < it4; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    trisp2[i + it2, j] = temp4[i, j];
                }
            }
            pointsplit = new double[q, m];
            for (int i = 0; i < q; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    pointsplit[i, j] = pointsplit1[i, j];
                }
            }
        }
        public void GetSquare()
        {
            double tempr, tempy;
            sqr = new double[q * 3, 3];
            leng = q * 3;
            int ity = 0, td;
            double[,] taboo = new double[q, 3];
            for (int k = 0; k < leng; k += 3)
            {
                tempr = l * 2;
                sqr[k, 0] = nxa;
                sqr[k, 1] = nya;
                sqr[k, 2] = nza;
                if (k >= 3)
                {
                    taboo[ity, 0] = sqr[k + 1, 0] = sqr[k - 1, 0];
                    taboo[ity, 1] = sqr[k + 1, 1] = sqr[k - 1, 1];
                    taboo[ity, 2] = sqr[k + 1, 2] = sqr[k - 1, 2];
                }
                else if (k < 2)
                {
                    taboo[ity, 0] = sqr[leng - 1, 0] = sqr[k + 1, 0] = pointsplit[0, 0];
                    taboo[ity, 1] = sqr[leng - 1, 1] = sqr[k + 1, 1] = pointsplit[0, 1];
                    taboo[ity, 2] = sqr[leng - 1, 2] = sqr[k + 1, 2] = pointsplit[0, 2];
                }
                ity++;
                for (int i = 0; i < q; i++)
                {
                    td = 1;
                    for (int e = 0; e < ity; e++) { if (pointsplit[i, 0] == taboo[e, 0] && pointsplit[i, 1] == taboo[e, 1] && pointsplit[i, 2] == taboo[e, 2]) td = 0; }
                    if (td == 1)
                    {
                        if (sqr[k + 1, 0] != pointsplit[i, 0] || sqr[k + 1, 1] != pointsplit[i, 1] || sqr[k + 1, 2] != pointsplit[i, 2])
                        {

                            tempy = Math.Sqrt(Math.Pow((sqr[k + 1, 0] - pointsplit[i, 0]), 2) + Math.Pow((sqr[k + 1, 1] - pointsplit[i, 1]), 2) + Math.Pow((sqr[k + 1, 2] - pointsplit[i, 2]), 2));

                            if (tempr >= tempy)
                            {

                                if (k >= 3)
                                {
                                    if (sqr[k - 2, 0] != pointsplit[i, 0] || sqr[k - 2, 1] != pointsplit[i, 1] || sqr[k - 2, 2] != pointsplit[i, 2])
                                    {
                                        tempr = tempy;
                                        sqr[k + 2, 0] = pointsplit[i, 0];
                                        sqr[k + 2, 1] = pointsplit[i, 1];
                                        sqr[k + 2, 2] = pointsplit[i, 2];
                                    }
                                }
                                else if (k < 3)
                                {
                                    tempr = tempy;
                                    sqr[k + 2, 0] = pointsplit[i, 0];
                                    sqr[k + 2, 1] = pointsplit[i, 1];
                                    sqr[k + 2, 2] = pointsplit[i, 2];
                                }
                            }
                        }
                    }
                }
            }          
        }
        public void WriteFileSplit(string path1, string path2)
        {
            int pp = 1, pn = 0;
            using (StreamWriter sw = new StreamWriter(path1))
            {
                sw.WriteLine("solid Cube");
                for (int i = 0; i < leng; i++)
                {
                    if (pp == 1)
                    {
                        sw.Write("facet normal");
                        sw.WriteLine();
                        sw.WriteLine("outer loop");
                    }
                    sw.Write("vertex");
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(Convert.ToString("\t" + sqr[i, j].ToString(CultureInfo.InvariantCulture)));
                    }
                    sw.Write("\n");
                    if (pp == 3)
                    {
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }
                    pp++;
                    if ((i + 1) % 3 == 0) pp = 1;
                    if ((i + 1) % 6 == 0) pn++;
                }
                for (int i = 0; i < itr1; i++)
                {
                    if (pp == 1)
                    {
                        sw.Write("facet normal");
                        sw.WriteLine();
                        sw.WriteLine("outer loop");
                    }
                    sw.Write("vertex");
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(Convert.ToString("\t" + trisp1[i, j].ToString(CultureInfo.InvariantCulture)));
                    }
                    sw.Write("\n");
                    if (pp == 3)
                    {
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }
                    pp++;
                    if ((i + 1) % 3 == 0) pp = 1;
                    if ((i + 1) % 6 == 0) pn++;
                }
                sw.WriteLine("endsolid Cube");
                sw.Close();
            }
            pp = 1;
            pn = 0;
            using (StreamWriter sw = new StreamWriter(path2))
            {
                sw.WriteLine("solid Cube");
                for (int i = 0; i < leng; i++)
                {
                    if (pp == 1)
                    {
                        sw.Write("facet normal");
                        sw.WriteLine();
                        sw.WriteLine("outer loop");
                    }
                    sw.Write("vertex");
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(Convert.ToString("\t" + sqr[i, j].ToString(CultureInfo.InvariantCulture)));
                    }
                    sw.Write("\n");
                    if (pp == 3)
                    {
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }
                    pp++;
                    if ((i + 1) % 3 == 0) pp = 1;
                    if ((i + 1) % 6 == 0) pn++;
                }
                for (int i = 0; i < itr2; i++)
                {
                    if (pp == 1)
                    {
                        sw.Write("facet normal");
                        sw.WriteLine();
                        sw.WriteLine("outer loop");
                    }
                    sw.Write("vertex");
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(Convert.ToString("\t" + trisp2[i, j].ToString(CultureInfo.InvariantCulture)));
                    }
                    sw.Write("\n");
                    if (pp == 3)
                    {
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }
                    pp++;
                    if ((i + 1) % 3 == 0) pp = 1;
                    if ((i + 1) % 6 == 0) pn++;
                }
                sw.WriteLine("endsolid Cube");
                sw.Close();
            }
            Console.WriteLine("Files succes created");
        }
    }
    public class Shere : IObgect
    {
        //1
        public double r, _r, a, a1, x, y, z, rad, Pi, b, b1;
        public int n1, m, n, g, horizont, pm;
        const int three = 3;
        public double[,] points;
        public double[,] triangles;
        //2
        public int todo = 1, q = 0, p = 0, itr1 = 0, itr2 = 0;
        public double[,] pointsplit1;
        public double[,] pointsplit;
        public double[,] pointsall;
        public double[,] trisp1;
        public double[,] trisp2;
        public double[,] sqr;
        public double nxa, nya, nza, nxb, nyb, nzb, A, B, C, D, T, x1, y1, z1, L = 0, M = 0, N = 0, mx, my, mz, l1, l2 = 0, l3 = 0;
        private int leng;
        public Shere(double _radius, double _x, double _y, double _z, int _m, int _horizont)
        {
            m = _m;
            horizont = _horizont;
            r = _radius;
            _r = -r;
            x = _x;
            y = _y;
            z = _z;
            g = horizont + 1;
            n = g + 1;
            b = 0;
            b1 = 360 / m;
            a = 0;
            a1 = 180 / g;
            Pi = 3.141592653589793238462643383279502884;
            rad = 180 / Pi;
            n1 = (int)(m * horizont + 2);
            pm = 1;
            points = new double[n, three];
            triangles = new double[n1, three];          
        }
        public void GetPoints()
        {
            for (int i = 0; i < n; i++)
            {
                points[i, 0] = (Math.Sin(a / rad) * r);
                points[i, 1] = 0;
                points[i, 2] = (Math.Cos(a / rad) * r);

                if (a == 180)
                {
                    points[i, 0] = 0;
                }
                if (a == 90)
                {
                    points[i, 2] = 0;
                }
                a += a1;
            }
        }

        public void GetTriangles()
        {
            triangles[0, 0] = 0;
            triangles[0, 1] = 0;
            triangles[0, 2] = points[0, 2];
            triangles[(int)(n1 - 1), 0] = 0;
            triangles[(int)(n1 - 1), 1] = 0;
            triangles[(int)(n1 - 1), 2] = points[(int)(n - 1), 2];

            for (int i = 1; i < g; i++)
            {
                double beta0 = 0;
                for (int k = pm; k < pm + m; k++)
                {
                    triangles[k, 0] = (Math.Cos(beta0 / rad) * points[i, 0]);
                    triangles[k, 1] = (Math.Sin(beta0 / rad) * points[i, 0]);
                    triangles[k, 2] = points[i, 2];
                    if (beta0 == 90 || beta0 == 270)
                    {
                        triangles[k, 0] = 0;
                    }
                    if (beta0 == 180 || beta0 == 360)
                    {
                        triangles[k, 1] = 0;
                    }
                    beta0 += b1;
                }
                pm += m;
            }

            for (int i = 0; i < n1; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    triangles[i, j] = Math.Round(triangles[i, j], 6);

                }
            }
        }

        public void WriteFile(string path)
        {
            using (StreamWriter sw = new StreamWriter(path))
            {
                sw.WriteLine("solid Sphere");
                for (int i = 0; i < n1; i++)
                {

                    if (i >= 1 && i < m + 1)
                    {
                        sw.Write("facet normal" + "\t");
                        sw.Write("\n");
                        sw.Write("outer loop" + "\n");
                        sw.Write("vertex");

                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[0, j].ToString(CultureInfo.InvariantCulture));

                        }

                        sw.Write("\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[i, j].ToString(CultureInfo.InvariantCulture));

                        }
                        sw.Write("\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            if (i < m)
                            {
                                sw.Write("\t" + triangles[i + 1, j].ToString(CultureInfo.InvariantCulture));
                            }
                            if (i == m)
                            {
                                sw.Write("\t" + triangles[1, j].ToString(CultureInfo.InvariantCulture));
                            }
                        }

                        sw.Write("\n");
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");

                    }

                    if (i > 0 && i < n1 - (m + 1))
                    {
                        sw.Write("facet normal" );
                        sw.Write("\n");
                        sw.Write("outer loop" + "\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[i, j].ToString(CultureInfo.InvariantCulture));
                        }
                        sw.Write("\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[i + 1, j].ToString(CultureInfo.InvariantCulture));
                        }
                        sw.Write("\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[i + m, j].ToString(CultureInfo.InvariantCulture));
                        }
                        sw.Write("\n");
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");



                        sw.Write("facet normal" );
                        sw.Write("\n");
                        sw.Write("outer loop" + "\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[i, j].ToString(CultureInfo.InvariantCulture));
                        }
                        sw.Write("\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[i + (m - 1), j].ToString(CultureInfo.InvariantCulture) );
                        }
                        sw.Write("\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[i + m, j].ToString(CultureInfo.InvariantCulture));
                        }
                        sw.Write("\n");
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }

                    if (i < n1 - 1 && i > n1 - (m + 2))
                    {
                        sw.Write("facet normal" );
                        sw.Write("\n");
                        sw.Write("outer loop" + "\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[n1 - 1, j].ToString(CultureInfo.InvariantCulture));
                        }
                        sw.Write("\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            sw.Write("\t" + triangles[i, j].ToString(CultureInfo.InvariantCulture));
                        }
                        sw.Write("\n");
                        sw.Write("vertex");
                        for (int j = 0; j < 3; j++)
                        {
                            if (i < n1 - 2)
                            {
                                sw.Write("\t" + triangles[i + 1, j].ToString(CultureInfo.InvariantCulture));
                            }
                            if (i == n1 - 2)
                            {
                                sw.Write("\t" + triangles[n1 - (m + 1), j].ToString(CultureInfo.InvariantCulture));
                            }
                        }

                        sw.Write("\n");
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }

                }
                sw.WriteLine("endsolid Sphere");
                sw.Close();
            }
            Console.WriteLine("File succes created");
        }
        public void ReadFile(string path)
        {
            int count = System.IO.File.ReadAllLines(path).Length;          
            count = ((count - 2) / 7) * 3;         
            n1 = count;
            triangles = new double[n1, three];
            pointsplit1 = new double[n1, three];


            using (StreamReader sr = new StreamReader(path))
            {
                string line, temps;
                int i = 0, tr = 0;
               

                while ((line = sr.ReadLine()) != null) 
                {


                    if (line.Contains('v'))
                    {
                        tr = 0;
                        char[] delimiterChars = { '\t' };

                        string text = line;

                        string[] words = text.Split(delimiterChars);
                        foreach (string word in words)
                        {
                            if (word != "vertex")
                            {                                
                                temps = word.Replace('.', ',');                               
                                triangles[i, tr] = Convert.ToDouble(temps);                               
                                tr++;
                            }
                        }
                        i++;                                              
                    }
                }
            }
            
        }
        public void GetSplit(double _nxa, double _nya, double _nza, double _nxb, double _nyb, double _nzb)
        {
            nxa = _nxa;
            nya = _nya;
            nza = _nza;

            nxb = _nxb;
            nyb = _nyb;
            nzb = _nzb;
            int it1 = 0, it2 = 0, it3 = 0, it4 = 0, chk = 0;
            double vx, vy, vz, nxc, nyc, nzc, length, lb1, lc1, lb2, lc2, b1, b2, b3, c1, c2, c3;

            double[,] temp1 = new double[n1, 3];
            double[,] temp2 = new double[n1, 3];
            double[,] temp3 = new double[n1, 3];
            double[,] temp4 = new double[n1, 3];

            vx = nxb - nxa;
            vy = nyb - nya;
            vz = nzb - nza;
            length = Math.Sqrt(Math.Pow((vx), 2) + Math.Pow((vy), 2) + Math.Pow((vz), 2));
            vx = vx / length;
            vy = vy / length;
            vz = vz / length;
            nxc = nxa - vx * length;
            nyc = nya - vy * length;
            nzc = nza - vz * length;

            pointsall = this.triangles;
            A = nxb - nxa;
            B = nyb - nya;
            C = nzb - nza;
            D = A * (-nxa) + B * (-nya) + C * (-nza);
            for (int i = 0; i < n1; i += 3)
            {
                chk = 0;
                //it1 = 0;
                //it2 = 0;

                b1 = Math.Sqrt(Math.Pow((nxb - pointsall[i, 0]), 2) + Math.Pow((nyb - pointsall[i, 1]), 2) + Math.Pow((nzb - pointsall[i, 2]), 2));
                c1 = Math.Sqrt(Math.Pow((nxc - pointsall[i, 0]), 2) + Math.Pow((nyc - pointsall[i, 1]), 2) + Math.Pow((nzc - pointsall[i, 2]), 2));

                b2 = Math.Sqrt(Math.Pow((nxb - pointsall[i + 1, 0]), 2) + Math.Pow((nyb - pointsall[i + 1, 1]), 2) + Math.Pow((nzb - pointsall[i + 1, 2]), 2));
                c2 = Math.Sqrt(Math.Pow((nxc - pointsall[i + 1, 0]), 2) + Math.Pow((nyc - pointsall[i + 1, 1]), 2) + Math.Pow((nzc - pointsall[i + 1, 2]), 2));

                b3 = Math.Sqrt(Math.Pow((nxb - pointsall[i + 2, 0]), 2) + Math.Pow((nyb - pointsall[i + 2, 1]), 2) + Math.Pow((nzb - pointsall[i + 2, 2]), 2));
                c3 = Math.Sqrt(Math.Pow((nxc - pointsall[i + 2, 0]), 2) + Math.Pow((nyc - pointsall[i + 2, 1]), 2) + Math.Pow((nzc - pointsall[i + 2, 2]), 2));
                if (b1 < c1 && b2 < c2 && b3 < c3)
                {
                    temp3[it3, 0] = pointsall[i, 0];
                    temp3[it3, 1] = pointsall[i, 1];
                    temp3[it3, 2] = pointsall[i, 2];
                    it3++;
                    temp3[it3, 0] = pointsall[i + 1, 0];
                    temp3[it3, 1] = pointsall[i + 1, 1];
                    temp3[it3, 2] = pointsall[i + 1, 2];
                    it3++;
                    temp3[it3, 0] = pointsall[i + 2, 0];
                    temp3[it3, 1] = pointsall[i + 2, 1];
                    temp3[it3, 2] = pointsall[i + 2, 2];
                    it3++;

                }
                else if (b1 > c1 && b2 > c2 && b3 > c3)
                {
                    temp4[it4, 0] = pointsall[i, 0];
                    temp4[it4, 1] = pointsall[i, 1];
                    temp4[it4, 2] = pointsall[i, 2];
                    it4++;
                    temp4[it4, 0] = pointsall[i + 1, 0];
                    temp4[it4, 1] = pointsall[i + 1, 1];
                    temp4[it4, 2] = pointsall[i + 1, 2];
                    it4++;
                    temp4[it4, 0] = pointsall[i + 2, 0];
                    temp4[it4, 1] = pointsall[i + 2, 1];
                    temp4[it4, 2] = pointsall[i + 2, 2];
                    it4++;
                }

                for (int j = 0; j < 2; j++)
                {
                    for (int k = 1; k < 3; k++)
                    {
                        if (j != k)
                        {
                            L = pointsall[i + j, 0] - pointsall[i + k, 0];
                            M = pointsall[i + j, 1] - pointsall[i + k, 1];
                            N = pointsall[i + j, 2] - pointsall[i + k, 2];
                            l1 = Math.Sqrt(Math.Pow(L, 2) + Math.Pow(M, 2) + Math.Pow(N, 2));
                            x1 = pointsall[i + j, 0];
                            y1 = pointsall[i + j, 1];
                            z1 = pointsall[i + j, 2];
                            if ((A * L + B * M + C * N) != 0)
                            {
                                T = -((A * x1 + B * y1 + C * z1 + D) / (A * L + B * M + C * N));
                                mx = L * T + x1;
                                my = M * T + y1;
                                mz = N * T + z1;
                                l2 = Math.Sqrt(Math.Pow((pointsall[i + j, 0] - mx), 2) + Math.Pow((pointsall[i + j, 1] - my), 2) + Math.Pow((pointsall[i + j, 2] - mz), 2));
                                l3 = Math.Sqrt(Math.Pow((pointsall[i + k, 0] - mx), 2) + Math.Pow((pointsall[i + k, 1] - my), 2) + Math.Pow((pointsall[i + k, 2] - mz), 2));
                                todo = 1;
                                for (int e = 0; e < q + 1; e++) { if (pointsplit1[e, 0] == mx && pointsplit1[e, 1] == my && pointsplit1[e, 2] == mz) todo = 0; }

                                bool bi = false;
                                double dep = l1 - (l2 + l3);
                                if (dep > 0 && dep < 0.001) { bi = true; }
                                else if (dep < 0 && dep > -0.001) { bi = true; }
                                else if (dep == 0 ) { bi = true; }

                                if (bi && todo == 1)
                                {
                                    pointsplit1[q, 0] = mx;
                                    pointsplit1[q, 1] = my;
                                    pointsplit1[q, 2] = mz;
                                    ++q;
                                }
                               
                                if (bi)
                                {
                                    lb1 = Math.Sqrt(Math.Pow((nxb - pointsall[i + j, 0]), 2) + Math.Pow((nyb - pointsall[i + j, 1]), 2) + Math.Pow((nzb - pointsall[i + j, 2]), 2));
                                    lc1 = Math.Sqrt(Math.Pow((nxc - pointsall[i + j, 0]), 2) + Math.Pow((nyc - pointsall[i + j, 1]), 2) + Math.Pow((nzc - pointsall[i + j, 2]), 2));

                                    lb2 = Math.Sqrt(Math.Pow((nxb - pointsall[i + k, 0]), 2) + Math.Pow((nyb - pointsall[i + k, 1]), 2) + Math.Pow((nzb - pointsall[i + k, 2]), 2));
                                    lc2 = Math.Sqrt(Math.Pow((nxc - pointsall[i + k, 0]), 2) + Math.Pow((nyc - pointsall[i + k, 1]), 2) + Math.Pow((nzc - pointsall[i + k, 2]), 2));

                                    if ((lb1 < lc1 && lb2 > lc2)|| (lb1 < lc1 && lb2 == lc2))
                                    {
                                        if (chk == 0)
                                        {
                                            //1
                                            temp1[it1, 0] = pointsall[i + j, 0];
                                            temp1[it1, 1] = pointsall[i + j, 1];
                                            temp1[it1, 2] = pointsall[i + j, 2];
                                            ++it1;
                                            temp1[it1, 0] = mx;
                                            temp1[it1, 1] = my;
                                            temp1[it1, 2] = mz;
                                            ++it1;
                                            //2
                                            temp2[it2, 0] = pointsall[i + k, 0];
                                            temp2[it2, 1] = pointsall[i + k, 1];
                                            temp2[it2, 2] = pointsall[i + k, 2];
                                            ++it2;
                                            temp2[it2, 0] = mx;
                                            temp2[it2, 1] = my;
                                            temp2[it2, 2] = mz;
                                            ++it2;
                                            chk++;
                                        }
                                        else if (chk == 1)
                                        {
                                            if (temp1[it1 - 2, 0] == pointsall[i + j, 0] && temp1[it1 - 2, 1] == pointsall[i + j, 1] && temp1[it1 - 2, 2] == pointsall[i + j, 2])
                                            {
                                                //1
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                //2
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                                temp2[it2, 0] = temp2[it2 - 3, 0];
                                                temp2[it2, 1] = temp2[it2 - 3, 1];
                                                temp2[it2, 2] = temp2[it2 - 3, 2];
                                                ++it2;
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                                temp2[it2, 0] = pointsall[i + k, 0];
                                                temp2[it2, 1] = pointsall[i + k, 1];
                                                temp2[it2, 2] = pointsall[i + k, 2];
                                                ++it2;
                                            }
                                            else
                                            {
                                                //1
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                temp1[it1, 0] = temp1[it1 - 3, 0];
                                                temp1[it1, 1] = temp1[it1 - 3, 1];
                                                temp1[it1, 2] = temp1[it1 - 3, 2];
                                                ++it1;
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                temp1[it1, 0] = pointsall[i + j, 0];
                                                temp1[it1, 1] = pointsall[i + j, 1];
                                                temp1[it1, 2] = pointsall[i + j, 2];
                                                ++it1;
                                                //2
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                            }
                                        }
                                    }
                                    else if ((lb1 > lc1 && lb2 < lc2)|| (lb1 == lc1 && lb2 < lc2))
                                    {
                                        if (chk == 0)
                                        {
                                            //1
                                            temp1[it1, 0] = pointsall[i + k, 0];
                                            temp1[it1, 1] = pointsall[i + k, 1];
                                            temp1[it1, 2] = pointsall[i + k, 2];
                                            ++it1;
                                            temp1[it1, 0] = mx;
                                            temp1[it1, 1] = my;
                                            temp1[it1, 2] = mz;
                                            ++it1;
                                            //2
                                            temp2[it2, 0] = pointsall[i + j, 0];
                                            temp2[it2, 1] = pointsall[i + j, 1];
                                            temp2[it2, 2] = pointsall[i + j, 2];
                                            ++it2;
                                            temp2[it2, 0] = mx;
                                            temp2[it2, 1] = my;
                                            temp2[it2, 2] = mz;
                                            ++it2;
                                            chk++;

                                        }
                                        else if (chk == 1)
                                        {
                                            if (temp1[it1 - 2, 0] == pointsall[i + k, 0] && temp1[it1 - 2, 1] == pointsall[i + k, 1] && temp1[it1 - 2, 2] == pointsall[i + k, 2])
                                            {
                                                //1
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                //2
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                                temp2[it2, 0] = temp2[it2 - 3, 0];
                                                temp2[it2, 1] = temp2[it2 - 3, 1];
                                                temp2[it2, 2] = temp2[it2 - 3, 2];
                                                ++it2;
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                                temp2[it2, 0] = pointsall[i + j, 0];
                                                temp2[it2, 1] = pointsall[i + j, 1];
                                                temp2[it2, 2] = pointsall[i + j, 2];
                                                ++it2;
                                            }
                                            else
                                            {
                                                //1
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                temp1[it1, 0] = temp1[it1 - 3, 0];
                                                temp1[it1, 1] = temp1[it1 - 3, 1];
                                                temp1[it1, 2] = temp1[it1 - 3, 2];
                                                ++it1;
                                                temp1[it1, 0] = mx;
                                                temp1[it1, 1] = my;
                                                temp1[it1, 2] = mz;
                                                ++it1;
                                                temp1[it1, 0] = pointsall[i + k, 0];
                                                temp1[it1, 1] = pointsall[i + k, 1];
                                                temp1[it1, 2] = pointsall[i + k, 2];
                                                ++it1;
                                                //2
                                                temp2[it2, 0] = mx;
                                                temp2[it2, 1] = my;
                                                temp2[it2, 2] = mz;
                                                ++it2;
                                            }

                                        }

                                    }

                                }

                            }
                        }
                    }
                }
            }
            itr1 = it1 + it3;
            itr2 = it2 + it4;
            trisp1 = new double[itr1, three];
            trisp2 = new double[itr2, three];
            for (int i = 0; i < it1; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    trisp1[i, j] = temp1[i, j];
                }
            }
            for (int i = 0; i < it3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    trisp1[i + it1, j] = temp3[i, j];
                }
            }
            for (int i = 0; i < it2; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    trisp2[i, j] = temp2[i, j];
                }
            }
            for (int i = 0; i < it4; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    trisp2[i+ it2, j] = temp4[i, j];
                }
            }
            pointsplit = new double[q, 3];
            for (int i = 0; i < q; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    pointsplit[i, j] = pointsplit1[i, j];
                }
            }
        }


        public void GetSquare()
        {
            double tempr, tempy;
            sqr = new double[q * 3, 3];
            leng = q * 3;
            int ity = 0,td;
            double[,] taboo = new double[q , 3];
            for (int k = 0; k < leng; k += 3)   
            {
                tempr = r * 2;
                sqr[k, 0] = nxa;
                sqr[k, 1] = nya;
                sqr[k, 2] = nza;
                if (k >= 3)
                {
                    taboo[ity, 0] = sqr[k + 1, 0] = sqr[k - 1, 0];
                    taboo[ity, 1] = sqr[k + 1, 1] = sqr[k - 1, 1];
                    taboo[ity, 2] = sqr[k + 1, 2] = sqr[k - 1, 2];                           
                }
                else if (k < 2)
                {
                    taboo[ity, 0] = sqr[leng - 1, 0] = sqr[k + 1, 0] = pointsplit[0, 0];
                    taboo[ity, 1] = sqr[leng - 1, 1] = sqr[k + 1, 1] = pointsplit[0, 1];
                    taboo[ity, 2] = sqr[leng - 1, 2] = sqr[k + 1, 2] = pointsplit[0, 2];                               
                }
                ity++;                
                for (int i = 0; i < q; i++)
                {
                    td = 1;
                    for (int e = 0; e < ity; e++) { if (pointsplit[i, 0] == taboo[e, 0] && pointsplit[i, 1] == taboo[e, 1] && pointsplit[i, 2] == taboo[e, 2]) td = 0; }
                    if (td == 1)
                    { 
                        if (sqr[k + 1, 0] != pointsplit[i, 0] || sqr[k + 1, 1] != pointsplit[i, 1] || sqr[k + 1, 2] != pointsplit[i, 2])
                        {

                            tempy = Math.Sqrt(Math.Pow((sqr[k + 1, 0] - pointsplit[i, 0]), 2) + Math.Pow((sqr[k + 1, 1] - pointsplit[i, 1]), 2) + Math.Pow((sqr[k + 1, 2] - pointsplit[i, 2]), 2));

                            if (tempr >= tempy)
                            {

                                if (k >= 3)
                                {
                                    if (sqr[k - 2, 0] != pointsplit[i, 0] || sqr[k - 2, 1] != pointsplit[i, 1] || sqr[k - 2, 2] != pointsplit[i, 2])
                                    {
                                        tempr = tempy;
                                        sqr[k + 2, 0] = pointsplit[i, 0];
                                        sqr[k + 2, 1] = pointsplit[i, 1];
                                        sqr[k + 2, 2] = pointsplit[i, 2];
                                    }
                                }
                                else if (k < 3)
                                {
                                    tempr = tempy;
                                    sqr[k + 2, 0] = pointsplit[i, 0];
                                    sqr[k + 2, 1] = pointsplit[i, 1];
                                    sqr[k + 2, 2] = pointsplit[i, 2];
                                }
                            }
                        }
                    }
                }
            }            
        }

        public void WriteFileSplit(string path1, string path2)
        {
            int pp = 1, pn = 0;
            using (StreamWriter sw = new StreamWriter(path1))
            {
                sw.WriteLine("solid Cube");
                for (int i = 0; i < leng; i++)
                {
                    if (pp == 1)
                    {
                        sw.Write("facet normal");
                        sw.WriteLine();
                        sw.WriteLine("outer loop");
                    }
                    sw.Write("vertex");
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(Convert.ToString("\t" + sqr[i, j].ToString(CultureInfo.InvariantCulture)));
                    }
                    sw.Write("\n");
                    if (pp == 3)
                    {
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }
                    pp++;
                    if ((i + 1) % 3 == 0) pp = 1;
                    if ((i + 1) % 6 == 0) pn++;
                }
                for (int i = 0; i < itr1; i++)
                {
                    if (pp == 1)
                    {
                        sw.Write("facet normal");
                        sw.WriteLine();
                        sw.WriteLine("outer loop");
                    }
                    sw.Write("vertex");
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(Convert.ToString("\t" + trisp1[i, j].ToString(CultureInfo.InvariantCulture)));
                    }
                    sw.Write("\n");
                    if (pp == 3)
                    {
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }
                    pp++;
                    if ((i + 1) % 3 == 0) pp = 1;
                    if ((i + 1) % 6 == 0) pn++;
                }
                sw.WriteLine("endsolid Cube");
                sw.Close();
            }
            pp = 1;
            pn = 0;
            using (StreamWriter sw = new StreamWriter(path2))
            {
                sw.WriteLine("solid Cube");
                for (int i = 0; i < leng; i++)
                {
                    if (pp == 1)
                    {
                        sw.Write("facet normal");
                        sw.WriteLine();
                        sw.WriteLine("outer loop");
                    }
                    sw.Write("vertex");
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(Convert.ToString("\t" + sqr[i, j].ToString(CultureInfo.InvariantCulture)));
                    }
                    sw.Write("\n");
                    if (pp == 3)
                    {
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }
                    pp++;
                    if ((i + 1) % 3 == 0) pp = 1;
                    if ((i + 1) % 6 == 0) pn++;
                }
                for (int i = 0; i < itr2; i++)
                {
                    if (pp == 1)
                    {
                        sw.Write("facet normal");
                        sw.WriteLine();
                        sw.WriteLine("outer loop");
                    }
                    sw.Write("vertex");
                    for (int j = 0; j < 3; j++)
                    {
                        sw.Write(Convert.ToString("\t" + trisp2[i, j].ToString(CultureInfo.InvariantCulture)));
                    }
                    sw.Write("\n");
                    if (pp == 3)
                    {
                        sw.WriteLine("endloop");
                        sw.WriteLine("endfacet");
                    }
                    pp++;
                    if ((i + 1) % 3 == 0) pp = 1;
                    if ((i + 1) % 6 == 0) pn++;
                }
                sw.WriteLine("endsolid Cube");
                sw.Close();
            }

            Console.WriteLine("Files succes created");
        }
    }
    class Program
    {
        static void Main(string[] args)
        {
            //Cube
            //1
            Cube cube = new Cube(10.0, 0.0, 0.0, 0.0);
            cube.GetPoints();
            cube.GetTriangles();
            cube.WriteFile(@"C:\\Users\\monst\\Desktop\\STl\\Cube.stl");
            //2
            cube.GetSplit(7, 7, 5, 19, 19, 5);
            cube.GetSquare();
            cube.WriteFileSplit(@"C:\\Users\\monst\\Desktop\\STl\\CubeSplit1.stl", @"C:\\Users\\monst\\Desktop\\STl\\CubeSplit2.stl");
            //Shere
            //1
            Shere shere = new Shere(100.0, 0.0, 0.0, 0.0, 24, 24);
            shere.GetPoints();
            shere.GetTriangles();
            shere.WriteFile(@"C:\\Users\\monst\\Desktop\\STl\\Sphere.stl");
            //2
            shere.ReadFile(@"C:\\Users\\monst\\Desktop\\STl\\Sphere.stl");
            //shere.GetSplit(60, 60, 10, 0, 0, 110);
            //shere.GetSplit(0, 0, 0, 10, 10, 110);
            //shere.GetSplit(10, 10, 0, 100, 10, 10);
            shere.GetSplit(10, 10, 80, 100, 10, 10);
            shere.GetSquare();
            shere.WriteFileSplit(@"C:\\Users\\monst\\Desktop\\STl\\SphereSplit1.stl", @"C:\\Users\\monst\\Desktop\\STl\\SphereSplit2.stl");
            ////
            Console.ReadKey();
        }
    }
}

using Microsoft.VisualBasic.Devices;
using System.Collections.Generic;
using System.Drawing.Drawing2D;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace Geodesic3D
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private double[,,] GenerateGeodesicVertices(int frequency, double radius)
        {
            // Calculate the golden ratio
            double phi = (1 + Math.Sqrt(5)) / 2;
            // Initialize the vertices array
            double[,,] vertices = new double[frequency + 1, frequency * 2 + 1, 3];
            // Generate the vertices
            for (int i = 0; i <= frequency; i++)
            {
                double y = i * radius / frequency - radius / 2;
                double phi1 = Math.Asin(y / radius);
                for (int j = 0; j <= frequency * 2; j++)
                {
                    double theta = j * Math.PI / frequency;
                    double phi2 = theta / phi;
                    double x = radius * Math.Cos(phi1) * Math.Cos(phi2);
                    double z = radius * Math.Cos(phi1) * Math.Sin(phi2);
                    vertices[i, j, 0] = x;
                    vertices[i, j, 1] = y;
                    vertices[i, j, 2] = z;
                }
            }
            // Generate the triangles
            List<int[]> triangles = new List<int[]>();
            for (int i = 0; i < frequency; i++)
            {
                for (int j = 0; j < frequency * 2; j++)
                {
                    int a = i * (frequency * 2 + 1) + j;
                    int b = a + 1;
                    int c = (i + 1) * (frequency * 2 + 1) + j;
                    triangles.Add(new int[] { a, b, c });
                    a = (i + 1) * (frequency * 2 + 1) + j;
                    b = b = i * (frequency * 2 + 1) + j + 1;
                    c = a + 1;
                    triangles.Add(new int[] { a, b, c });
                }
            }
            // Flatten the vertices and triangles into arrays
            double[] flattenedVertices = new double[(frequency + 1) * (frequency * 2 + 1) * 3];
            int[] flattenedTriangles = new int[triangles.Count * 3];
            int vertexIndex = 0;
            foreach (double d in vertices)
            {
                flattenedVertices[vertexIndex++] = d;
            }
            int triangleIndex = 0;
            foreach (int[] t in triangles)
            {
                flattenedTriangles[triangleIndex++] = t[0];
                flattenedTriangles[triangleIndex++] = t[1];
                flattenedTriangles[triangleIndex++] = t[2];
            }
            return vertices;
        }

        private double[,] CalculateTransform(int width, int height)
        {
            double aspectRatio = (double)width / height;
            double fieldOfView = 60.0;
            double near = 0.1;
            double far = 100.0;

            // Compute the perspective projection matrix
            double[,] projection = new double[4, 4];
            double fovRadians = Math.PI * fieldOfView / 180.0;
            double yScale = 1.0 / Math.Tan(fovRadians / 2.0);
            double xScale = yScale / aspectRatio;
            double zRange = far - near;
            projection[0, 0] = xScale;
            projection[1, 1] = yScale;
            projection[2, 2] = -(far + near) / zRange;
            projection[2, 3] = -2.0 * far * near / zRange;
            projection[3, 2] = -1.0;

            // Compute the view matrix
            double[,] view = new double[4, 4];
            view[0, 0] = 1.0;
            view[1, 1] = 1.0;
            view[2, 2] = 1.0;
            view[3, 3] = 1.0;
            view[3, 2] = -10.0;

            // Compute the world matrix
            double[,] world = new double[4, 4];
            world[0, 0] = 1.0;
            world[1, 1] = 1.0;
            world[2, 2] = 1.0;
            world[3, 3] = 1.0;

            // Calculate the transform matrix
            double[] flattenedworld = world.Cast<double>().ToArray();
            double[] flattenedview = view.Cast<double>().ToArray();
            double[] views = MatrixMultiply(projection, flattenedview);
            double[,] viewt = Reshape(views, 4, 4);
            double[] trans = MatrixMultiply(viewt, flattenedworld);
            double[,] transform = Reshape(trans, 4, 4);

            return transform;

        }


        private void RenderGeodesic()
        {
            // Create a bitmap to draw the geodesic on
            Bitmap bitmap = new Bitmap(pictureBox1.Width, pictureBox1.Height);

            // Get the graphics object to draw on the bitmap
            Graphics graphics = Graphics.FromImage(bitmap);

            // Set the background color of the bitmap to white
            graphics.Clear(Color.White);
            int width = pictureBox1.Width;
            int height = pictureBox1.Height;
            // Call GenerateGeodesicVertices() to generate the triangles
            double[,,] vertices = GenerateGeodesicVertices(10, 5);
            // Define the triangles array
            int[,] triangles = {
                { 0, 11, 5 },
                { 0, 5, 1 },
                { 0, 1, 7 },
                { 0, 7, 10 },
                { 0, 10, 11 },
                { 1, 5, 9 },
                { 5, 11, 4 },
                { 11, 10, 2 },
                { 10, 7, 6 },
                { 7, 1, 8 },
                { 3, 9, 4 },
                { 3, 4, 2 },
                { 3, 2, 6 },
                { 3, 6, 8 },
                { 3, 8, 9 },
                { 4, 9, 5 },
                { 2, 4, 11 },
                { 6, 2, 10 },
                { 8, 6, 7 },
                { 9, 8, 1 }
            };
            // Calculate the transformation matrix to transform the geodesic vertices to screen coordinates
            double[,] transform = CalculateTransform(width, height);


            // Draw each triangle in the geodesic
            for (int i = 0; i < triangles.GetLength(0); i++)
            {
                double ax = vertices[triangles[i, 0], triangles[i, 1], 0];
                double ay = vertices[triangles[i, 0], triangles[i, 1], 1];
                double az = vertices[triangles[i, 0], triangles[i, 1], 2];
                double bx = vertices[triangles[i, 2], triangles[i, 3], 0];
                double by = vertices[triangles[i, 2], triangles[i, 3], 1];
                double bz = vertices[triangles[i, 2], triangles[i, 3], 2];
                double cx = vertices[triangles[i, 4], triangles[i, 5], 0];
                double cy = vertices[triangles[i, 4], triangles[i, 5], 1];
                double cz = vertices[triangles[i, 4], triangles[i, 5], 2];
                DrawTriangle(pictureBox1.CreateGraphics(), transform, ax, ay, az, bx, by, bz, cx, cy, cz);
            }
            // Set the picture box's image to the bitmap
            pictureBox1.Image = bitmap;
        }

        private double[] MatrixMultiply(double[,] matrix, double[] vector)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            //if (vector.Length != cols)
            //    throw new ArgumentException("Vector length must match matrix column count.");

            double[] result = new double[rows];
            for (int i = 0; i < rows; i++)
            {
                double sum = 0;
                for (int j = 0; j < cols; j++)
                {
                    sum += matrix[i, j] * vector[j];
                }
                result[i] = sum;
            }

            return result;

        }

        private void DrawTriangle(Graphics graphics, double[,] transform, double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
        {
            // Apply the transform to the vertices
            double[] a = MatrixMultiply(transform, new double[] { ax, ay, az, 1 });
            double[] b = MatrixMultiply(transform, new double[] { bx, by, bz, 1 });
            double[] c = MatrixMultiply(transform, new double[] { cx, cy, cz, 1 });

            // Check if the triangle is facing the camera
            if (CalculateNormal(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2])[2] > 0)
            {
                // Convert the vertices to screen coordinates
                double sx = a[0] / a[3];
                double sy = a[1] / a[3];
                double sz = a[2] / a[3];
                double bx2 = b[0] / b[3];
                double by2 = b[1] / b[3];
                double cx2 = c[0] / c[3];
                double cy2 = c[1] / c[3];

                // Scale and translate the screen coordinates
                double width = pictureBox1.Width / 2;
                double height = pictureBox1.Height / 2;
                sx = sx * width + width;
                sy = -sy * height + height;
                bx2 = bx2 * width + width;
                by2 = -by2 * height + height;
                cx2 = cx2 * width + width;
                cy2 = -cy2 * height + height;

                // Draw the triangle
                graphics.DrawLine(Pens.Black, (float)sx, (float)sy, (float)bx2, (float)by2);
                graphics.DrawLine(Pens.Black, (float)sx, (float)sy, (float)cx2, (float)cy2);
                graphics.DrawLine(Pens.Black, (float)bx2, (float)by2, (float)cx2, (float)cy2);
            }
        }

        private double[] CalculateNormal(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz)
        {
            // Calculate the vectors between the vertices
            double[] ab = new double[] { bx - ax, by - ay, bz - az };
            double[] ac = new double[] { cx - ax, cy - ay, cz - az };

            // Calculate the normal vector
            double nx = ab[1] * ac[2] - ab[2] * ac[1];
            double ny = ab[2] * ac[0] - ab[0] * ac[2];
            double nz = ab[0] * ac[1] - ab[1] * ac[0];

            // Normalize the normal vector
            double length = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            return new double[] { nx / length, ny / length, nz / length };
        }

        private double[,] Reshape(double[] array, int rows, int columns)
        {
            //if (array.Length != rows * columns)
            //    throw new ArgumentException("Array length must be equal to rows * columns.");

            double[,] result = new double[rows, columns];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    result[i, j] = array[i * columns + j];
                }
            }
            return result;
        }

        //private void DrawGoedesic(int frequency, double radius)
        //{
        //    // Generate the vertices
        //    double[,,] vertices = GenerateGeodesicVertices(frequency, radius);
        //    // Create a new bitmap for the picture box
        //    Bitmap bitmap = new Bitmap(pictureBox1.Width, pictureBox1.Height);
        //    // Get a graphics object for the bitmap
        //    Graphics graphics = Graphics.FromImage(bitmap);
        //    // Clear the bitmap
        //    graphics.Clear(Color.White);
        //    // Set the projection matrix
        //    double fov = Math.PI / 4;
        //    double aspectRatio = (double)pictureBox1.Width / (double)pictureBox1.Height;
        //    double near = 0.1;
        //    double far = 1000;
        //    double[,] projection = new double[,]
        //    {
        //        { 1 / (Math.Tan(fov / 2) * aspectRatio), 0, 0, 0 },
        //        { 0, 1 / Math.Tan(fov / 2), 0, 0 },
        //        { 0, 0, (far + near) / (far - near), 1 },
        //        { 0, 0, -(2 * far * near) / (far - near), 0 }
        //    };
        //    // Set the view matrix
        //    double[,] view = new double[,]
        //    {
        //        { 1, 0, 0, 0 },
        //        { 0, 1, 0, 0 },
        //        { 0, 0, 1, -radius * 2 },
        //        { 0, 0, 0, 1 }
        //    };
        //    // Set the world matrix
        //    double[,] world = new double[,]
        //    {
        //        { 1, 0, 0, 0 },
        //        { 0, 1, 0, 0 },
        //        { 0, 0, 1, 0 },
        //        { 0, 0, 0, 1 }
        //    };
        //    // Calculate the transform matrix
        //    double[] flattenedworld = world.Cast<double>().ToArray();
        //    double[] flattenedview = view.Cast<double>().ToArray();
        //    double[] views = MatrixMultiply(projection, flattenedview);
        //    double[,] viewt = Reshape(views, 4, 4);
        //    double[] trans = MatrixMultiply(viewt, flattenedworld);
        //    double[,] transform = Reshape(trans, 4, 4);

        //    // Draw the triangles
        //    for (int i = 0; i < frequency; i++)
        //    {
        //        for (int j = 0; j < frequency * 2; j++)
        //        {
        //            int a = i * (frequency * 2 + 1) + j;
        //            int b = a + 1;
        //            int c = (i + 1) * (frequency * 2 + 1) + j;
        //            DrawTriangle(graphics, transform, vertices[a, 0], vertices[a, 1], vertices[a, 2], vertices[b, 0], vertices[b, 1], vertices[b, 2], vertices[c, 0], vertices[c, 1], vertices[c, 2]);
        //            a = (i + 1) * (frequency * 2 + 1) + j;
        //            b = b = i * (frequency * 2 + 1) + j + 1;
        //            c = a + 1;
        //            DrawTriangle(graphics, transform, vertices[a, 0], vertices[a, 1], vertices[a, 2], vertices[b, 0], vertices[b, 1], vertices[b, 2], vertices[c, 0], vertices[c, 1], vertices[c, 2]);
        //        }
        //    }
        //    // Set the picture box's image to the bitmap
        //    pictureBox1.Image = bitmap;
        //}

        private void Form1_Load(object sender, EventArgs e)
        {
            // Initialize the 3D array to store the coordinates of the vertices
            double[,,] vertices = GenerateGeodesicVertices(10, 5);
            // Define the radius of the geodesic
            //double radius = 180;
            // Define the number of subdivisions
            //int subdivisions = 20;
            // Define the rotation angle in radians
            double angle = Math.PI / 4;
            // Define the rotation axis as a unit vector
            double[] axis = new double[] { 0, 1, 0 };
            double axisLength = Math.Sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
            axis[0] /= axisLength;
            axis[1] /= axisLength;
            axis[2] /= axisLength;
            // Define the rotation matrix
            double[,] rotationMatrix = new double[3, 3];
            rotationMatrix[0, 0] = Math.Cos(angle) + axis[0] * axis[0] * (1 - Math.Cos(angle));
            rotationMatrix[0, 1] = axis[0] * axis[1] * (1 - Math.Cos(angle)) - axis[2] * Math.Sin(angle);
            rotationMatrix[0, 2] = axis[0] * axis[2] * (1 - Math.Cos(angle)) + axis[1] * Math.Sin(angle);
            rotationMatrix[1, 0] = axis[1] * axis[0] * (1 - Math.Cos(angle)) + axis[2] * Math.Sin(angle);
            rotationMatrix[1, 1] = Math.Cos(angle) + axis[1] * axis[1] * (1 - Math.Cos(angle));
            rotationMatrix[1, 2] = axis[1] * axis[2] * (1 - Math.Cos(angle)) - axis[0] * Math.Sin(angle);
            rotationMatrix[2, 0] = axis[2] * axis[0] * (1 - Math.Cos(angle)) - axis[1] * Math.Sin(angle);
            rotationMatrix[2, 1] = axis[2] * axis[1] * (1 - Math.Cos(angle)) + axis[0] * Math.Sin(angle);
            rotationMatrix[2, 2] = Math.Cos(angle) + axis[2] * axis[2] * (1 - Math.Cos(angle));
            // Loop through each vertex
            for (int i = 0; i < vertices.GetLength(0); i++)
            {
                for (int j = 0; j < vertices.GetLength(1); j++)
                {
                    // Get the coordinates of the vertex
                    double x = vertices[i, j, 0];
                    double y = vertices[i, j, 1];
                    double z = vertices[i, j, 2];

                    // Apply the rotation transformation
                    double newX = rotationMatrix[0, 0] * x + rotationMatrix[0, 1] * y + rotationMatrix[0, 2] * z;
                    double newY = rotationMatrix[1, 0] * x + rotationMatrix[1, 1] * y + rotationMatrix[1, 2] * z;
                    double newZ = rotationMatrix[2, 0] * x + rotationMatrix[2, 1] * y + rotationMatrix[2, 2] * z;

                    // Update the coordinates of the vertex
                    vertices[i, j, 0] = newX;
                    vertices[i, j, 1] = newY;
                    vertices[i, j, 2] = newZ;

                }
            }

        }

        private void pictureBox1_Paint(object sender, PaintEventArgs e)
        {
            // Render the geodesic when the picture box is painted
            RenderGeodesic();
        }
    }
}
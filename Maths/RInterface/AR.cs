#region Imporar

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

#endregion

namespace Maths {
    
    internal class AR {

        #region Fields

        private double[] arr;

        #endregion

        #region Constructor

        internal AR(int length) {
            this.arr = new double[length+1];
        }

        internal AR(double[] arr) {
            this.arr = new double[arr.Length + 1];
            for (int i = 0; i < arr.Length; i++) { this.arr[i+1] = arr[i]; }
        }

        internal AR(int[] arr) {
            this.arr = new double[arr.Length + 1];
            for (int i = 0; i < arr.Length; i++) { this.arr[i+1] = (double)arr[i]; }
        }

        #endregion

        #region Properties
        
        internal int Length {
            get { return arr.Length - 1; }
        }

        internal double this[int index] {
            get {
                if (index < 1 || index > Length) { throw new Exception("Error. Out of bounds"); }
                return this.arr[index]; 
            }
            set {
                if (index < 1 || index > Length) { throw new Exception("Error. Out of bounds"); }
                this.arr[index] = value;
            }
        }

        #endregion

        #region Overrides

        public override string ToString() {
            string ar = "";
            for (int i = 1; i <= Length; i++) { ar = ar + this[i] + " "; }
            return ar + "\n";
        }

        #endregion

        #region Operators

        public static AR operator +(AR ar, double val) {
            AR arSum = new AR(ar.Length);
            for (int i = 1; i <= ar.Length; i++) { arSum[i] = ar[i] + val; }
            return arSum;
        }

        public static AR operator -(AR ar, double val) {
            AR arSub = new AR(ar.Length);
            for (int i = 1; i <= ar.Length; i++) { arSub[i] = ar[i] - val; }
            return arSub;
        }

        public static AR operator *(AR ar, double val) {
            AR arPro = new AR(ar.Length);
            for (int i = 1; i <= ar.Length; i++) { arPro[i] = ar[i] * val; }
            return arPro;
        }

        public static AR operator /(AR ar, double val) {
            if (val == 0) { throw new Exception("Error. Division by zero"); }
            AR arDiv = new AR(ar.Length);
            for (int i = 1; i <= ar.Length; i++) { arDiv[i] = ar[i] / val; }
            return arDiv;
        }

        public static AR operator +(AR ar1, AR ar2) {
            if(ar1.Length != ar2.Length) { throw new Exception("Error. Different size"); }
            AR arSum = new AR(ar1.Length);
            for (int i = 1; i <= ar1.Length; i++) { arSum[i] = ar1[i] + ar2[i]; }
            return arSum;    
        }

        public static AR operator -(AR ar1, AR ar2) {
            if (ar1.Length != ar2.Length) { throw new Exception("Error. Different size"); }
            AR arSub = new AR(ar1.Length);
            for (int i = 1; i <= ar1.Length; i++) { arSub[i] = ar1[i] - ar2[i]; }
            return arSub;
        }

        public static AR operator *(AR ar1, AR ar2) {
            if (ar1.Length != ar2.Length) { throw new Exception("Error. Different size"); }
            AR arPro = new AR(ar1.Length);
            for (int i = 1; i <= ar1.Length; i++) { arPro[i] = ar1[i] * ar2[i]; }
            return arPro;
        }

        public static AR operator /(AR ar1, AR ar2) {
            if (ar1.Length != ar2.Length) { throw new Exception("Error. Different size"); }
            AR arDiv = new AR(ar1.Length);
            for (int i = 1; i <= ar1.Length; i++) { arDiv[i] = ar1[i] / ar2[i]; }
            return arDiv;
        }

        #endregion

        #region Static Methods

        internal static AR Abs(AR ar) {
            AR arAbs = ar.Clone();
            for (int i = 1; i <= ar.Length; i++) { arAbs[i] = Math.Abs(ar[i]); }
            return arAbs;
        }
        
        internal static AR C(AR ar1, AR ar2) {
            AR conc = new AR(ar1.Length + ar2.Length);
            Copy(ar1, conc, 1, ar1.Length, 1, ar1.Length);
            Copy(ar2, conc, 1, ar2.Length, ar1.Length + 1, conc.Length);
            return conc;
        }

        internal static AR C(AR ar, double x) {
            AR conc = new AR(ar.Length + 1);
            Copy(ar, conc, 1, ar.Length, 1, ar.Length);
            conc[conc.Length] = x;
            return conc;
        }

        internal static double Sum(AR ar1, AR ar2) {
            if (ar1.Length != ar2.Length) { throw new Exception("Error. Different size"); }
            double sum = 0;
            for (int i = 1; i <= ar1.Length; i++) { sum = sum + ar1[i] + ar2[i]; }
            return sum;
        }

        internal static double Sum(AR ar) {
            return ar.Sum();
        }
        
        internal static void Copy(AR arFrom, AR arTo, int iniFrom, int endFrom, int iniTo, int endTo) {
            if (endFrom < iniFrom || endTo < iniTo) { throw new Exception("Error. end cannot be less than ini"); }
            if (endFrom - iniFrom != endTo - iniTo) { throw new Exception("Error. intervals from and to are different"); }
            int diff = iniTo - iniFrom;
            for (int i = iniTo; i <= endTo; i++) { arTo[i] = arFrom[i - diff]; }
        }

        internal static AR Rep(double val, int n) {
            AR rep = new AR(n);
            for (int i = 1; i <= rep.Length; i++) { rep[i] = val; }
            return rep;
        }
        
        #endregion

        #region Internal Methods

        internal double[] ToArray() {
            double[] zeroIndexedArray = new double[Length];
            for (int i = 1; i <= Length; i++) { zeroIndexedArray[i - 1] = this[i];  }
            return zeroIndexedArray;
        }

        internal List<double> ToList() {
            List<double> list = new List<double>();
            for (int i = 1; i <= Length; i++) { list.Add(this[i]); }
            return list;
        }
        
        internal double[] GetArr() { return this.arr; }
        
        internal AR Clone() {
            AR arClone = new AR(this.Length);
            AR.Copy(this, arClone, 1, this.Length, 1, this.Length);
            return arClone;
        } 

        internal AR SubAr(int ini, int end) {
            int diff = end - ini;
            AR subar = new AR(diff + 1);
            for (int i = ini; i <= end; i++) { subar[i-ini+1] = this[i]; }
            return subar;
        }

        internal double Sum(int ini, int end) {
            double sum = 0;
            for (int i = ini; i <= end; i++) { sum += this[i]; }
            return sum;
        }

        internal double Sum() {
            return Sum(1, this.Length);
        }

        internal double Min() {
            double min = double.MaxValue;
            for (int i = 1; i <= Length; i++) { if (this[i] < min) min = this[i]; }
            return min;
        }

        internal double Max() {
            double max = double.MinValue;
            for (int i = 1; i <= Length; i++) { if (this[i] > max) max = this[i]; }
            return max;
        }
        
        internal void ElimNegatives() {
            for (int i = 1; i <= this.Length; i++) { if (this[i] < 0) this[i] = 0; }
        }

        #endregion

        #region Statistics

        internal double Mean() {
            if (Length == 0) { return 0.0; }
            return this.Sum() / this.Length;
        }

        internal double Variance() {
            int n = this.Length;
            if (n == 1) { return 0; }
            double Sum = 0.0, SumSquares = 0.0;

            for (int i = 1; i <= n; i++) {
                Sum += this[i];
                SumSquares += this[i] * this[i];
            }
            return (n * SumSquares - (Sum * Sum)) / (n * n - 1);
        }

        internal double Sd() {
            return Math.Sqrt(Variance());
        }

        #endregion

    }
}

#region Imports

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

#endregion

namespace Maths {
    
    internal class DF {

        #region Fields

        private List<List<double>> cols;

        #endregion

        #region Constructors

        internal DF() {
            this.cols = new List<List<double>>();
        }

        internal DF(int nCol) {
            this.cols = new List<List<double>>();
            for (int i = 0; i < nCol; i++) { this.cols.Add(new List<double>()); }
        }


        internal DF(double[][] mat) : this() {
            int c = mat.Length;
            int r = mat[0].Length;
            for (int i = 0; i < c; i++) { this.cols.Add(mat[i].ToList<double>()); }
        }

        internal DF(double[] arr, bool isRow) : this() {
            if (isRow) {
                for (int i = 0; i < arr.Length; i++) {
                    this.cols.Add(new List<double>());
                    this.cols[i].Add(arr[i]);
                }
            }
            else { 
                this.cols.Add(arr.ToList<double>()); 
            }
        }

        internal DF(AR ts, bool isRow) : this(ts.ToArray(), isRow) {
        }

        internal DF(DF df) : this() {
            List<List<double>> dfCols = df.GetCols();
            this.cols = new List<List<double>>();
            for (int i = 0; i < dfCols.Count; i++) { 
                this.cols.Add(new List<double>());
                for (int j = 0; j < dfCols[i].Count; j++) { this.cols[i].Add(dfCols[i][j]); }
            }
        }
        
        #endregion

        #region Properties

        internal List<List<double>> GetCols() {
            return cols; 
        }
        
        internal int NCol {
            get { return cols.Count; }
        }

        internal int NRow {
            get {
                if (cols.Count == 0) { return 0; }
                return cols[0].Count; 
            }
        }

        internal double this[int row, int col] {
            get {
                if (col < 1 || col > NCol) { throw new Exception("Col out of bounds"); }
                if (row < 1 || row > NRow) { throw new Exception("Row out of bounds"); }
                return cols[col-1][row-1]; 
            }
            set {
                if (col < 1 || col > NCol) { throw new Exception("Col out of bounds"); }
                if (row < 1 || row > NRow) { throw new Exception("Row out of bounds"); }
                cols[col-1][row-1] = value;
            }
        }
        
        #endregion

        #region Overrides

        public override string ToString() {
            string df = "";
            for (int i = 1; i <= NRow; i++) {
                for (int j = 1; j <= NCol; j++) {
                    df = df + this[i, j] + "\t";
                }
                df = df + "\n";
            }
            return df;
        }

        internal DF Clone() { 
            DF dfClone = new DF(this);
            return dfClone;
        }

        public List<double> GetColumn(int index) {
            return this.cols[index-1];
        }

        public List<double> GetRow(int index) {
            List<double> row = new List<double>();
            for (int j = 1; j <= NCol; j++) {
                row.Add(this[index, j]);
            }
            return row;
        }
       
        #endregion

        #region Operators

        public static DF operator +(DF df, double val) {
            DF dfRes = df.Clone();
            for (int i = 1; i <= df.NRow; i++) {
                for (int j = 1; j <= df.NCol; j++) {
                    dfRes[i,j] = df[i,j] + val;
                }   
            }
            return dfRes;
        }

        public static DF operator -(DF df, double val) {
            DF dfRes = df.Clone();
            for (int i = 1; i <= df.NRow; i++) {
                for (int j = 1; j <= df.NCol; j++) {
                    dfRes[i, j] = df[i, j] - val;
                }
            }
            return dfRes;
        }

        public static DF operator *(DF df, double val) {
            DF dfRes = df.Clone();
            for (int i = 1; i <= df.NRow; i++) {
                for (int j = 1; j <= df.NCol; j++) {
                    dfRes[i, j] = df[i, j] * val;
                }
            }
            return dfRes;
        }

        public static DF operator /(DF df, double val) {
            DF dfRes = df.Clone();
            for (int i = 1; i <= df.NRow; i++) {
                for (int j = 1; j <= df.NCol; j++) {
                    dfRes[i, j] = df[i, j] / val;
                }
            }
            return dfRes;
        }

        public static DF operator +(DF df1, DF df2) {
            if (df1.NRow != df2.NRow || df1.NCol != df2.NCol) { throw new Exception("Error. Different dimensions"); }
            DF dfRes = df1.Clone();
            for (int i=1;i<=dfRes.NRow;i++) {
                for (int j = 1; j <= dfRes.NCol; j++) {
                    dfRes[i, j] += df2[i, j]; 
                }
            }
            return dfRes;
        }

        public static DF operator -(DF df1, DF df2) {
            if (df1.NRow != df2.NRow || df1.NCol != df2.NCol) { throw new Exception("Error. Different dimensions"); }
            DF dfRes = df1.Clone();
            for (int i = 1; i <= dfRes.NRow; i++) {
                for (int j = 1; j <= dfRes.NCol; j++) {
                    dfRes[i, j] -= df2[i, j];
                }
            }
            return dfRes;
        }

        public static DF operator *(DF df1, DF df2) {
            if (df1.NRow != df2.NRow || df1.NCol != df2.NCol) { throw new Exception("Error. Different dimensions"); }
            DF dfRes = df1.Clone();
            for (int i = 1; i <= dfRes.NRow; i++) {
                for (int j = 1; j <= dfRes.NCol; j++) {
                    dfRes[i, j] *= df2[i, j];
                }
            }
            return dfRes;
        }


        public static DF operator /(DF df1, DF df2) {
            if (df1.NRow != df2.NRow || df1.NCol != df2.NCol) { throw new Exception("Error. Different dimensions"); }
            DF dfRes = df1.Clone();
            for (int i = 1; i <= dfRes.NRow; i++) {
                for (int j = 1; j <= dfRes.NCol; j++) {
                    dfRes[i, j] /= df2[i, j];
                }
            }
            return dfRes;
        }

        #endregion
        
        #region Static Methods

        internal static DF AddColumn(DF df, List<double> col) {
            DF dfRes = df.Clone();
            dfRes.AddColumn(col);
            return dfRes;
        }

        internal static DF AddColumn(DF df, AR ts) {
            List<double> col = ts.ToList();
            DF dfRes = df.Clone();
            dfRes.AddColumn(col);
            return dfRes;
        }
   
        internal static DF AddRow(DF df1, List<double> row) {
            DF dfRes = df1.Clone();
            dfRes.AddRow(row);
            return dfRes;
        }

        internal static DF AddRow(DF df, AR ts) {
            List<double> row = ts.ToList();
            DF dfRes = AddRow(df, row);
            return dfRes;
        }

        internal static DF CBind(DF df1, DF df2) {
            DF dfRes = df1.Clone();
            for (int i = 1; i <= df2.NCol; i++) { dfRes.AddColumn(df2.GetColumn(i)); }
            return dfRes;
        }
        
        internal static DF RBind(DF df1, DF df2) {
            DF dfBind = df1.Clone();
            for (int i = 1; i <= df2.NRow; i++) {
                List<double> row = df2.GetRow(i);
                dfBind.AddRow(row);
            }
            return dfBind;
        }

        #endregion
        
        #region Internal Methods

        internal void AddColumn(List<double> col) {
            this.cols.Add(col);
        }

        internal void AddColumn(AR ts) {
            List<double> col = ts.ToList();
            this.AddColumn(col);
        }

        public void AddRow(List<double> row) {
            if (row.Count != NCol) { throw new Exception("Error. Wrong dimensions"); }
            for (int j = 1; j <= NCol; j++) {
                this.cols[j-1].Add(row[j-1]);
            }
        }

        internal void AddRow(AR ts) {
            List<double> row = ts.ToList();
            this.AddRow(row);
        }

        internal void CBind(DF df) {
            for (int i = 1; i <= df.NCol; i++) { this.AddColumn(df.GetColumn(i)); }
        }

        internal void RBind(DF df) {
            for (int i = 1; i <= df.NRow; i++) {
                List<double> row = df.GetRow(i);
                this.AddRow(row);
            }
        }

        internal DF SubDfCols(int ini, int end) {
            DF dfRes = new DF();
            for (int i = ini; i <= end; i++) {
                dfRes.AddColumn(this.GetColumn(i));
            }
            return dfRes;
        }

        internal DF SubDfCols(int[] sel) {
            DF dfRes = new DF();
            for (int i = 0; i < sel.Length; i++) {
                dfRes.AddColumn(this.GetColumn(sel[i]));
            }
            return dfRes;
        }
        
        internal DF SubDfRows(int ini, int end) {
            DF dfRes = new DF(this.NCol);
            for (int i = ini; i <= end; i++) {
                dfRes.AddRow(this.GetRow(i));
            }
            return dfRes;
        }

        internal DF SubDfRows(int[] sel) {
            DF dfRes = new DF(this.NCol);
            for (int i = 0; i < sel.Length; i++) {
                dfRes.AddRow(this.GetRow(sel[i]));
            }
            return dfRes;
        }

        internal double[][] ToColsArray() {
            double[][] colsArr = new double[NCol][];
            for(int i=0;i<NCol;i++) {
                colsArr[i] = cols[i].ToArray();
            }
            return colsArr;
        }

        internal double[][] ToRowsArray() {
            double[][] rowsArr = new double[NRow][];
            for (int i = 0; i < NRow; i++) {
                rowsArr[i] = new double[NCol];
                for (int j = 0; j < NCol; j++) {
                    rowsArr[i][j] = cols[j][i];
                }
            }
            return rowsArr;
        }
        
        internal List<List<double>> ToList() {
            return cols;
        }
        
        #endregion

    }
}

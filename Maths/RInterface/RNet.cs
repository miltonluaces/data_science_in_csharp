#region Imports

using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using RDotNet;
using Maths;
using System.IO;

#endregion

namespace Maths {

    public class RNet {

        #region Fields

        private string path;
        private string version;
        private REngine engine;
        private SymbolicExpression ret;

        #endregion

        #region Constructor

        internal RNet(string path, string version) {
            this.path = path;
            this.version = version;

            REnv.GetInstance(path, version);
            this.engine = REngine.GetInstance();
            this.engine.AutoPrint = false;
            ret = engine.Evaluate("path='" + path + "'");
            engine.Evaluate("setwd('" + path + "')");

            this.engine = REngine.GetInstance();
            this.engine.AutoPrint = false;
            engine.Evaluate("path='" + path + "'");
            engine.Evaluate("setwd('" + path + "')");
        }

        private void SetupPath(string path, string rVersion) {
            var oldPath = System.Environment.GetEnvironmentVariable("PATH");
            var rPath = System.Environment.Is64BitProcess ? path + "R/" + rVersion + "/bin/x64" : path + "R/" + rVersion + "/bin/i386";
            if (!Directory.Exists(rPath)) throw new DirectoryNotFoundException(string.Format(" R.dll not found in : {0}", rPath));
            var newPath = string.Format("{0}{1}{2}", rPath, System.IO.Path.PathSeparator, oldPath);
            System.Environment.SetEnvironmentVariable("PATH", newPath);
        }

        #endregion

        #region Internal Methods

        internal void LoadLibrary(string libName) {
            engine.Evaluate("library(" + libName + ")");
        }

        internal void RequireLibrary(string libraryName) {
            CharacterVector package = engine.CreateCharacter(libraryName);
            engine.SetSymbol("package", package);
            engine.Evaluate("if(libName %in% rownames(installed.packages()))  do.call('library', list(package))  else { install.packages(package) do.call('library', list(package))");
        }

        internal void Execute(string callStr) {
            ret = engine.Evaluate(callStr);
        }

        internal void Execute(string callStr, IList<Params> pars) {
            SymbolicExpression se = null;
            for (int i = 0; i < pars.Count; i++) {
                switch (pars[i].type) {
                    case ParamType.Numeric:
                        se = engine.CreateNumeric(Convert.ToDouble(pars[i].value));
                        break;
                    case ParamType.NumericVector:
                        double[] arr = ((AR)(pars[i].value)).ToArray();
                        se = engine.CreateNumericVector(((IList)arr).Cast<double>());
                        break;
                    case ParamType.IntVector:
                        se = engine.CreateIntegerVector(((IList)pars[i].value).Cast<int>());
                        break;
                    case ParamType.String:
                        se = engine.CreateCharacter(pars[i].value.ToString());
                        break;
                    case ParamType.Bool:
                        bool[] bools = { (bool)pars[i].value };
                        se = engine.CreateLogicalVector(((IList)bools).Cast<bool>());
                        break;
                    case ParamType.BoolVector:
                        se = engine.CreateLogicalVector(((IList)pars[i].value).Cast<bool>());
                        break;
                    case ParamType.Dataframe:
                        double[][] rowsArr = ((DF)(pars[i].value)).ToRowsArray();
                        se = engine.CreateDataFrame((IEnumerable[])rowsArr);
                        break;
                    case ParamType.REnvironment:
                        se = engine.CreateEnvironment((REnvironment)pars[i].value);
                        break;
                    case ParamType.SymbolicExpression:
                        se = (SymbolicExpression)pars[i].value;
                        break;
                }
                engine.SetSymbol(pars[i].name, se);
            }
            ret = engine.Evaluate(callStr);
        }

        internal S4Object GetS4Object() {
            return ret.AsS4();
        }

        internal int[] GetIntReturn() {
            return ret.AsInteger().ToArray();
        }

        internal double[] GetDoubleReturn() {
            return ret.AsNumeric().ToArray();
        }

        internal string[] GetStringReturn() {
            return ret.AsCharacter().ToArray();
        }

        internal bool[] GetBoolReturn() {
            return ret.AsLogical().ToArray();
        }

        internal double[][] GetDataframeReturn() {
            DataFrame df = ret.AsDataFrame();
            double[][] dfArr = new double[df.ColumnCount][];
            for (int i = 0; i < df.ColumnCount; i++) {
                dfArr[i] = new double[df.RowCount];
                for (int j = 0; j < df.RowCount; j++) { dfArr[i][j] = Convert.ToDouble(df[i][j]); }
            }
            return dfArr;
        }

        internal SymbolicExpression GetReturn() {
            return ret;
        }

        internal REnvironment GetEnvironmentReturn() {
            REnvironment env = ret.AsEnvironment();
            return env;
        }

        internal double[][] ToMatrix(DataFrame df) {
            double[][] dfArr = new double[df.ColumnCount][];
            for (int i = 0; i < df.ColumnCount; i++) {
                dfArr[i] = new double[df.RowCount];
                for (int j = 0; j < df.RowCount; j++) { dfArr[i][j] = Convert.ToDouble(df[i][j]); }
            }
            return dfArr;
        }

        internal GenericVector GetGenericVectorReturn() {
            return ret.AsList();
        }

        internal void CreateR6Object(string className, string objName) {
            Execute(objName + " = " + className + "$new()");
            REnvironment env = GetEnvironmentReturn();
            engine.SetSymbol(objName, env);
        }

        internal void LoadObject(string objName) {
            REnvironment env = GetEnvironmentReturn();
            engine.SetSymbol(objName, env);
        }

        internal Dictionary<string, object> GetReturnObject(GenericVector gv) {
            Dictionary<string, object> res = new Dictionary<string, object>();
            for (int i = 0; i < gv.Length; i++) {
                switch (gv[i].Type) {
                    case RDotNet.Internals.SymbolicExpressionType.CharacterVector:
                        res.Add(gv.Names[i], gv[i].AsCharacter().ToArray());
                        break;
                    case RDotNet.Internals.SymbolicExpressionType.List:
                        res.Add(gv.Names[i], ToMatrix(gv[0].AsDataFrame()));
                        break;
                    case RDotNet.Internals.SymbolicExpressionType.NumericVector:
                        res.Add(gv.Names[i], gv[i].AsNumeric().ToArray());
                        break;
                    case RDotNet.Internals.SymbolicExpressionType.IntegerVector:
                        res.Add(gv.Names[i], gv[i].AsInteger().ToArray());
                        break;
                }
            }
            return res;
        }

        #endregion

        #region Diagnostics

        internal bool JITEnabled() {
            ret = engine.Evaluate("Sys.getenv('R_COMPILE_PKGS')");
            string rcp = GetStringReturn()[0];
            ret = engine.Evaluate("Sys.getenv('R_ENABLE_JIT')");
            string rej = GetStringReturn()[0];
            return (rcp == "TRUE" && rej == "3");
        }

        internal void SetRSeed(int seed) {
            IList<RNet.Params> pars = new List<RNet.Params>(); 
            pars.Add(new RNet.Params("seed", RNet.ParamType.Numeric,  seed));
            if (seed >= 0) { Execute("set.seed(seed)", pars); } else { Execute("set.seed(null)"); }
        }

        internal string GetWorkingDirectory() {
            ret = engine.Evaluate("getwd()");
            string wd = GetStringReturn()[0];
            return wd.ToString();
        }

        internal string GetPackageVersions() {
            ret = engine.Evaluate("packageVersion('MiscFunctions')");
            string mfv = GetStringReturn()[0];
            ret = engine.Evaluate("packageVersion('Multivariate')");
            string mvv = GetStringReturn()[0];
            return "MiscFunctions: " + mfv + " - Multivariate: " + mvv;
        }

        #endregion

        #region Private Methods

        private System.Collections.IEnumerable GetIEnumerable(List<List<double>> columns) {
            IEnumerable[] ieColumns = new IEnumerable[columns.Count];
            for (int i = 0; i < columns.Count; i++) {
                ieColumns[i] = columns[i].ToArray();
            }
            return ieColumns;
        }


        #endregion

        #region Enums

        internal enum ParamType { Numeric, NumericVector, IntVector, String, Bool, BoolVector, Dataframe, S4Object, REnvironment, SymbolicExpression };

        #endregion

        #region Inner Classes

        internal class Params {
            internal Params(string name, ParamType type, object value) {
                this.name = name;
                this.type = type;
                this.value = value;
            }
            public string name;
            public ParamType type;
            public object value;
        }

        #endregion

    }
}
#region Imports

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

#endregion

namespace Maths {

    public class REnv {

        #region Fields

        private static object sync = new object();
        private static REnv instance;

        #endregion

        #region Constructor

        /// <summary> Constructor, private due to the singleton </summary>
        private REnv() {

        }

        #endregion

        #region Singleton method

        /// <summary> Obtains the singleton instance </summary>
        public static REnv GetInstance(string path, string version) {
            if (REnv.instance == null) {
                lock (sync) {
                    if (REnv.instance == null) {
                        instance = new REnv();
                        instance.SetEnvironmentVariables(path, version);
                    }
                }
            }
            return REnv.instance;
        }

        #endregion

        #region Private Methods

        private void SetEnvironmentVariables(string path, string version) {
            var oldPath = System.Environment.GetEnvironmentVariable("PATH");
            var rPath = path + "R/R-" + version + "/bin/i386";
            var newPath = string.Format("{0}{1}{2}", rPath, System.IO.Path.PathSeparator, oldPath);
            System.Environment.SetEnvironmentVariable("PATH", newPath);
            string rHome = System.Environment.GetEnvironmentVariable("R_HOME");
            if (string.IsNullOrEmpty(rHome)) { rHome = path + "R/R-" + version; }
            System.Environment.SetEnvironmentVariable("R_HOME", rHome);
        }

        #endregion
    }
}
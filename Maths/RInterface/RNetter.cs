#region Imports

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

#endregion

namespace Maths {

    public interface RNetter {

        bool JITEnabled();
        string GetWorkingDirectory();
        void SetSeed(int seed);
        string GetPackageVersions();
    }
}

// Copyright (C) 2014 Alexey Matveichev
// Copyright (C) 2011-2013 OpenFOAM Foundation
//
// DISCLAIMER
// This source code is not approved or endorsed by OpenCFD Limited, producer and
// distributor of the OpenFOAM software via www.openfoam.com, and owner of the
// OPENFOAM(R)  and OpenCFD(R)  trade marks.
//
// ACKNOWLEDGEMENT
// OPENFOAM(R)  is a registered trade mark of OpenCFD Limited, producer and
// distributor of the OpenFOAM software via www.openfoam.com.
//
// DESCRIPTION
// Solves a transport equation for a passive scalar with diffusivity as a sum
// molecular of turbulent part. Turbulent part is calculated from turbulent
// viscosity and turbulent Schmidt number.

#include "fvCFD.H"
#include "fvIOoptionList.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
              - fvm::laplacian(DTt, T)
             ==
                fvOptions(T)
            );
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

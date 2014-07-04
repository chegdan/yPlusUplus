/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    yPlus

Description
    Calculates and reports yPlus for all wall patches, for the specified times.

    Default behaviour assumes operating in incompressible mode.
    Use the -compressible option for compressible RAS cases.

    Copyright 2013 Daniel Wei

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fluidThermo.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "wallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    argList::addBoolOption
    (
        "compressible",
        "calculate compressible y+"
    );

    argList::addBoolOption
    (
        "writeY",
        "write wall distance files"
    );

    argList::addBoolOption
    (
        "writeYPlus",
        "write yPlus files"
    );

    argList::addBoolOption
    (
        "uTauInfo",
        "Show uTau infor, used in channelFlow post-process"
    );

#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    const bool compressible = args.optionFound("compressible");
    const bool writeY = args.optionFound("writeY");
    const bool writeYPlus = args.optionFound("writeYPlus");
    const bool uTauInfo = args.optionFound("uTauInfo");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        //-Wall distance
        //-y has values on all mesh cells
        if (writeY)
        {
            fvMesh::readUpdateState state = mesh.readUpdate();
            if (timeI == 0 || state != fvMesh::UNCHANGED)
            {
                Info<< "Calculating wall distance\n" << endl;
                wallDist y(mesh, true);
                Info<< "Writing wall distance to field " << y.name() << nl << endl;
                y.write();
            }
        }

        //-wallDistanceOnPatch only has values on wall patch cells
        // volScalarField wallDistanceOnPatch
        // (
        //     IOobject
        //     (
        //         "wallDistanceOnPatch",
        //         runTime.timeName(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::NO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar("wallDistanceOnPatch", dimLength, 0.0)
        // );

        //-yPlus only has values on wall patch cells
        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        //-uTau only has values on wall patch cells
        volScalarField uTau
        (
            IOobject
            (
                "uTau",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("uTau", dimVelocity, 0.0)
        );

        Info<< "Reading field U\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volScalarField muEff
        (
            IOobject
            (
                "muEff",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("muEff", dimensionSet(1,-1,-1,0,0,0,0), 0.0)
        );
        volScalarField muLam
        (
            IOobject
            (
                "muLam",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("muLam", dimensionSet(1,-1,-1,0,0,0,0), 0.0)
        );
        volScalarField nuEff
        (
            IOobject
            (
                "nuEff",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("nuEff", dimensionSet(0,2,-1,0,0,0,0), 0.0)
        );
        volScalarField nuLam
        (
            IOobject
            (
                "nuLam",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("nuLam", dimensionSet(0,2,-1,0,0,0,0), 0.0)
        );

        if (compressible)
        {
            autoPtr<fluidThermo> pThermo
            (
                fluidThermo::New(mesh)
            );
            fluidThermo& thermo = pThermo();

            volScalarField rho
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                thermo.rho()
            );

            #include "compressibleCreatePhi.H"

            autoPtr<compressible::turbulenceModel> comModel
            (
                compressible::turbulenceModel::New
                (
                    rho,
                    U,
                    phi,
                    thermo
                )
            );

            muEff=comModel->muEff();
            muLam=comModel->mu();
        }
        else
        {
#           include "createPhi.H"

            singlePhaseTransportModel laminarTransport(U, phi);

            autoPtr<incompressible::turbulenceModel> incModel
            (
                incompressible::turbulenceModel::New(U, phi, laminarTransport)
            );

            nuEff=incModel->nuEff();
            nuLam=incModel->nu();
        }

        volScalarField::GeometricBoundaryField d = nearWallDist(mesh).y();

        const fvPatchList& patches = mesh.boundary();
        
        dimensionedScalar uTauAvg("uTauAvg", dimVelocity, 0);
        scalar nPatch = 0;

        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (isA<wallFvPatch>(currPatch))
            {
                if (compressible)
                {
                    autoPtr<basicThermo> pThermo
                    (
                        basicThermo::New(mesh)
                    );
                    basicThermo& thermo = pThermo();

                    volScalarField rho
                    (
                        IOobject
                        (
                            "rho",
                            runTime.timeName(),
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE
                        ),
                        thermo.rho()
                    );

                    // wallDistanceOnPatch.boundaryField()[patchi] = d[patchi];
                    yPlus.boundaryField()[patchi] =
                        d[patchi]
                       *sqrt
                        (
                            rho.boundaryField()[patchi]
                           *muEff.boundaryField()[patchi]
                           *mag(U.boundaryField()[patchi].snGrad())
                        )
                       /muLam.boundaryField()[patchi];

                    uTau.boundaryField()[patchi] =
                        sqrt
                        (
                            rho.boundaryField()[patchi]
                           *muEff.boundaryField()[patchi]
                           *mag(U.boundaryField()[patchi].snGrad())
                        );
                }
                else
                {
                    // wallDistanceOnPatch.boundaryField()[patchi] = d[patchi];
                    yPlus.boundaryField()[patchi] =
                        d[patchi]
                       *sqrt
                        (
                            nuEff.boundaryField()[patchi]
                           *mag(U.boundaryField()[patchi].snGrad())
                        )
                       /nuLam.boundaryField()[patchi];

                    uTau.boundaryField()[patchi] =
                        sqrt
                        (
                            nuEff.boundaryField()[patchi]
                           *mag(U.boundaryField()[patchi].snGrad())
                        );
                }

                const scalarField& Yp = yPlus.boundaryField()[patchi];

                // Mean uTau estimation 
                uTauAvg.value() += average(uTau.boundaryField()[patchi]);
                nPatch ++;

                Info<< "Patch " << patchi
                    << " named " << currPatch.name() << nl
                    << " d   : min: " << min(d[patchi])
                    << " max: " << max(d[patchi])
                    << " average: " << average(d[patchi]) << nl 
                    << " uTau: min: " << min(uTau.boundaryField()[patchi])
                    << " max: " << max(uTau.boundaryField()[patchi])
                    << " average: " << average(uTau.boundaryField()[patchi])<<nl 
                    << " y+  : min: " << min(Yp) << " max: " << max(Yp)
                    << " average: " << average(Yp) << nl << endl;

            }
        }

        // Average over all walls, used in channelFlow post-processing
        uTauAvg /= nPatch; 
        if (uTauInfo)
        {
            Info << "Average friction velocity uTau is: "
                 << uTauAvg.value()
                 << " (averaged over " << nPatch << " wall(s))"
                 << nl <<endl;
        }

        if (writeYPlus)
        {
            Info<< "Writing yPlus to field "
                << yPlus.name() << nl << endl;
            yPlus.write();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

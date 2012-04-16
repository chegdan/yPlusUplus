/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    plusPostRANS

Description

    calculates y+ and u+ fields for wall-bounded flows computed with
    one of the available low-Re RANS (no wall function!) turbulence 
    models. More specifically it

    :: 	calculates and outputs y+ (avg., min., max.) based on the 
	velocity gradient at the wall  

    ::	calculates and outputs the wall averaged friction velocity 

    ::  writes fields of y+ and U+ to the corresponding time directory

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "nearWallDist.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        fvMesh::readUpdateState state = mesh.readUpdate();

        wallDist y(mesh, true);

        if (timeI == 0 || state != fvMesh::UNCHANGED)
        {
            Info<< "Calculating wall distance\n" <<endl;
            Info<< "Writing wall distance to field " << y.name() << nl << endl;
            y.write();
        }

	#include "createFields.H"

        const fvPatchList& patches = mesh.boundary();
                     
	Info<< "Summary: " << nl << endl;

        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (isA<wallFvPatch>(currPatch))
            {
                yPlusTemp.boundaryField()[patchi] =
                    d[patchi]
                   *sqrt
                    (
                        RASModel->nu()().boundaryField()[patchi]
                       *mag(U.boundaryField()[patchi].snGrad())
                    )
                   /RASModel->nu()().boundaryField()[patchi];
                
		const scalarField& YpTemp = yPlusTemp.boundaryField()[patchi];

                uTau.boundaryField()[patchi] =
                    sqrt
                    (
		    RASModel->nu()()
                   *mag(U.boundaryField()[patchi].snGrad())
                    );

                Info<< "  y+ for Patch " << patchi
                    << " named " << currPatch.name() << ":" 
                    << " min: " << min(YpTemp) << " max: " << max(YpTemp)
		    << nl << endl;
            }
        }


const volVectorField& centers = mesh.C();
const surfaceVectorField& faceCenters = mesh.Cf();

	forAll(uTau, cellI){


		forAll(patches, patchi){
            		const fvPatch& currPatch = patches[patchi];
	
			label nFaces = mesh.boundaryMesh()[patchi].size();
        
			if(isA<wallFvPatch>(currPatch)){
				for(int facei = 0; facei<nFaces; facei++){

				//dimensionedScalar cellFaceDist("cellFaceDist",dimensionSet(0,0,0,0,0,0,0),scalar(1));
				scalar cellFaceDist ;

				cellFaceDist = Foam::sqrt(sqr(centers[cellI].x()-faceCenters.boundaryField()[patchi][facei].x()) + sqr(centers[cellI].y()-faceCenters.boundaryField()[patchi][facei].y())+ sqr(centers[cellI].z()-faceCenters.boundaryField()[patchi][facei].z()));

				//convert the y value for comparison
				scalar yTemp = y[cellI];

				scalar diffDist = abs(cellFaceDist - yTemp)/max(abs(cellFaceDist),SMALL);

				//compare the values
				//if( cellFaceDist == yTemp){ uTau[cellI] = uTau.boundaryField()[patchi][facei];	}
				if( diffDist <= 0.001){ uTau[cellI] = uTau.boundaryField()[patchi][facei];	}
	
					}
			}	
		}


	}
//uTau.write();
        yPlus = y.y() * uTau / RASModel->nu();

	dimensionedScalar UUnits ( "UUnits", dimensionSet(0,1,-1,0,0), 1.0 ); 

	uPlus = U / stabilise(uTau, SMALL*UUnits);
	
	/*forAll(uPlus, cellI){
	
        	uPlus[cellI].x() = U[cellI].x() / max(uTau[cellI],SMALL);
       		uPlus[cellI].y() = U[cellI].y() / max(uTau[cellI],SMALL);
        	uPlus[cellI].z() = U[cellI].z() / max(uTau[cellI],SMALL);

	}*/
       
        Info << "Writing yPlus and uPlus to corresponding fields." << nl <<endl;
        yPlus.write();
        uPlus.write();
	uTau.write();

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

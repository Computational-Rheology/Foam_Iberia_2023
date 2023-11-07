/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "newBoundaryConditionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newBoundaryConditionFvPatchScalarField::newBoundaryConditionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF) // Initialize parent class
{}


Foam::newBoundaryConditionFvPatchScalarField::newBoundaryConditionFvPatchScalarField
(
    const newBoundaryConditionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper) // Initialize parent class
{}


Foam::newBoundaryConditionFvPatchScalarField::newBoundaryConditionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false) // Initialize parent class
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())                    // If `value` is defined, initialize with what is defined in value
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(scalarField(p.size(), 1.23)); // Else, initialize the field with 1.23
    }
}


Foam::newBoundaryConditionFvPatchScalarField::newBoundaryConditionFvPatchScalarField
(
    const newBoundaryConditionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf) // Initialize parent class
{}


Foam::newBoundaryConditionFvPatchScalarField::newBoundaryConditionFvPatchScalarField
(
    const newBoundaryConditionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)     // Initialize parent class
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::newBoundaryConditionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
        
    // view values in watch for scalarField: *scalarField.v_@nPoints
    // view values in watch for vectorField: vectorField.v_.v_@nPoints
    
    // grab reference to fvPatch
    const fvPatch& p = patch(); 
    
    // Grab information from fvPatch
    
    const word& name = p.name();   // patch name

    const vectorField& Cf = p.Cf();  // (x,y,z) position of face center

    const Foam::vectorField& Sf = p.Sf(); // Face area vector

    const Foam::scalarField& magSf = p.magSf(); // Area of face
    
    vectorField nf = p.nf();    // face unit normal

    vectorField delta = p.delta(); // distance to the face projected in the normal direction
                                  // nHat*(nHat & (Cf() - Cn()))

    const scalarField& deltaCoefs = p.deltaCoeffs();    // reciprocal of face distance projected on the normal direction

    vectorField Cn = p.Cn();    // (x,y,z) of cell center of cell adjacent to the face

    scalarField cellCenterValues = patchInternalField(); // values at the cell center adjacent to patch

    const scalarField& domainValues = internalField(); // values of internalField

    scalarField normalGrad = snGrad(); // Surface normal gradient without non-orthogonal correction 
                                   // patch().deltaCoeffs()*(*this - patchInternalField());

    const scalarField& faceValues = *this;      // Face values

    // Look up patch of another field 

    const auto& U = db().lookupObject<volVectorField>("U");     // volVectorfield
 
    const auto& phi = db().lookupObject<surfaceScalarField>("phi"); // surfaceScalarField

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>("U");  // Same patch, different field

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    operator==
    (
        scalarField(p.size(), 20)       // fixedValue 10
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField with given weights
Foam::tmp<Foam::Field<Foam::scalar>>  Foam::newBoundaryConditionFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>& t
) const
{
    Info << "valueInternalCoeffs" << endl;
    return fixedValueFvPatchScalarField::valueInternalCoeffs(t);
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the value of this patchField with given weights
Foam::tmp<Foam::Field<Foam::scalar>> Foam::newBoundaryConditionFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>& t
) const
{
    Info << "valueBoundaryCoeffs" << endl;
    return fixedValueFvPatchScalarField::valueBoundaryCoeffs(t);
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
Foam::tmp<Foam::Field<Foam::scalar>> Foam::newBoundaryConditionFvPatchScalarField::gradientInternalCoeffs() const
{
    Info << "gradientInternalCoeffs" << endl;
    return fixedValueFvPatchScalarField::gradientInternalCoeffs();
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the gradient of this patchField
Foam::tmp<Foam::Field<Foam::scalar>> Foam::newBoundaryConditionFvPatchScalarField::gradientBoundaryCoeffs() const
{
    Info << "gradientBoundaryCoeffs" << endl;
    return fixedValueFvPatchScalarField::gradientBoundaryCoeffs();
}


void Foam::newBoundaryConditionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        newBoundaryConditionFvPatchScalarField
    );
}

// ************************************************************************* //

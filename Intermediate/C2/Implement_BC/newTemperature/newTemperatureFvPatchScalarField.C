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

#include "newTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newTemperatureFvPatchScalarField::newTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    TSwitch_(450),
    TCooling_(300),
    THeating_(600)

{}


Foam::newTemperatureFvPatchScalarField::newTemperatureFvPatchScalarField
(
    const newTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    TSwitch_(ptf.TSwitch_),
    TCooling_(ptf.TCooling_),
    THeating_(ptf.THeating_)
{}


Foam::newTemperatureFvPatchScalarField::newTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    TSwitch_(readScalar(dict.lookup("TSwitch"))),
    TCooling_(readScalar(dict.lookup("TCooling"))),
    THeating_(readScalar(dict.lookup("THeating")))

{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    

}


Foam::newTemperatureFvPatchScalarField::newTemperatureFvPatchScalarField
(
    const newTemperatureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    TSwitch_(tppsf.TSwitch_),
    TCooling_(tppsf.TCooling_),
    THeating_(tppsf.THeating_)
  
    
{}


Foam::newTemperatureFvPatchScalarField::newTemperatureFvPatchScalarField
(
    const newTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    TSwitch_(tppsf.TSwitch_),
    TCooling_(tppsf.TCooling_),
    THeating_(tppsf.THeating_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::newTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::newTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::newTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Acess the dimensioned internalField
    const DimensionedField<scalar, volMesh>& tInt = 
            this->internalField();

    //Calculate cell volume weighted average of the field
    const dimensionedScalar tAverage =
        tInt.weightedAverage(patch().boundaryMesh().mesh().V());

    //set the fixed value boundary condition according specified switch value
    //and calculated average value

    if(tAverage.value() < TSwitch_)
    {
        operator==(THeating_);
    }
    else
    {
        operator==(TCooling_); 
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::newTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("TSwitch") << TSwitch_ << token::END_STATEMENT << nl;
    os.writeKeyword("TCooling") << TCooling_ << token::END_STATEMENT << nl;
    os.writeKeyword("THeating") << THeating_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        newTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //

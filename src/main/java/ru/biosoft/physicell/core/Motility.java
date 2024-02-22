package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.List;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;

/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
public class Motility implements Cloneable
{
    public boolean isMotile;
    public double persistenceTime; // mean time to keep going in one direction before resampling for a new direction. 
    public double migrationSpeed; // migration speed along chosen direction, in absence of all other adhesive / repulsive forces 
    public double[] migrationBiasDirection = new double[0];; // a unit vector random motility is biased in this direction (e.g., chemotaxis)
    public double migrationBias; // how biased is motility if 0, completely random. if 1, deterministic along the bias vector 
    public boolean restrictTo2D; // if true, set random motility to 2D only. 
    public double[] motilityVector = new double[0];;
    public int chemotaxisIndex;
    public int chemotaxisDirection;
    public double[] chemotacticSensitivities = new double[0]; // advanced chemotaxis
    private String[] substrates;

    public Motility()
    {
        isMotile = false;
        persistenceTime = 1.0;
        migrationSpeed = 1.0;
        migrationBiasDirection = new double[3];
        migrationBias = 0.0;
        restrictTo2D = false;
        motilityVector = new double[3];
        chemotaxisIndex = -1;
        chemotaxisDirection = 1;
        substrates = new String[0];
    }

    void sync(Microenvironment m)
    {
        chemotacticSensitivities = VectorUtil.resize( chemotacticSensitivities, m.numberDensities(), 0 );
        substrates = m.densityNames;
    }

    public void setChemotacticSensitivity(int index, double value)
    {
        chemotacticSensitivities[index] = value;
    }

    public void disable()
    {
        isMotile = false;
        motilityVector = new double[3];
    }

    @Override
    public Motility clone()
    {
        try
        {
            Motility result = (Motility)super.clone();
            result.chemotacticSensitivities = this.chemotacticSensitivities.clone();
            result.motilityVector = this.motilityVector.clone();
            result.migrationBiasDirection = this.migrationBiasDirection.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }

    public String display()
    {
        if( !isMotile )
            return "Motility Disabled.\n--------------------------------";
        StringBuilder sb = new StringBuilder();
        sb.append( "Motility\n--------------------------------" );

        sb.append( "\n\tIn " + ( restrictTo2D ? "2D" : "3D" ) );
        sb.append( ", speed: " + migrationSpeed + " micron/min, bias: " + migrationBias + ", persistence: " + persistenceTime + " min" );
        if( chemotaxisIndex != -1 )
            sb.append( "\n\tChemotaxis along " + chemotaxisDirection + " * ( " + substrates[chemotaxisIndex] + " )" );

        List<String> chemoSenses = new ArrayList<>();
        for( int i = 0; i < chemotacticSensitivities.length; i++ )
        {
            double sens = chemotacticSensitivities[i];
            if( sens != 0 )
                chemoSenses.add( sens + " * ( " + substrates[i] + " )" );
        }

        if( !chemoSenses.isEmpty() )
        {
            sb.append( "\n\tChemotaxis along  " + chemoSenses.get( 0 ) );
            for( int i = 1; i < chemoSenses.size(); i++ )
                sb.append( "\n\t                 +" + chemoSenses.get( i ) );
        }
        return sb.toString();
    }
}
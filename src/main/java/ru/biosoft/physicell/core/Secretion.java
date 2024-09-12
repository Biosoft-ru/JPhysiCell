package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.List;

import ru.biosoft.physicell.biofvm.BasicAgent;
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
public class Secretion implements Cloneable
{
    public double[] secretionRates = new double[0];
    public double[] uptakeRates = new double[0];
    public double[] saturationDensities = new double[0];
    public double[] netExportRates = new double[0];
    private String[] substrates = new String[0];

    public void sync(Microenvironment m)
    {
        int size = m.numberDensities();
        substrates = m.densityNames;
        secretionRates = VectorUtil.resize( secretionRates, size );
        uptakeRates = VectorUtil.resize( uptakeRates, size );
        saturationDensities = VectorUtil.resize( saturationDensities, size, 1.0 );
        netExportRates = VectorUtil.resize( netExportRates, size );
    }

    public void advance(BasicAgent cell, Phenotype phenotype, double dt)
    {
        if( cell == null ) // if this phenotype is not associated with a cell, exit 
            return;

        // make sure the associated cell has the correct rate vectors 
        if( cell.secretionRates != secretionRates )
        {
            cell.secretionRates = secretionRates;//TODO: remove cell.secretionRate
            cell.uptakeRates = uptakeRates;
            cell.saturationDensities = saturationDensities;
            cell.netExportRates = netExportRates;
            cell.setTotalVolume( phenotype.volume.total );
            cell.setUptakeConstants( dt );
        }
        cell.simulateSecretionUptake( cell.getMicroenvironment(), dt );
    }

    public void setSecretionToZero()
    {
        for( int i = 0; i < secretionRates.length; i++ )
        {
            secretionRates[i] = 0.0;
            netExportRates[i] = 0.0;
        }
    }

    public void setUptakeToZero()
    {
        for( int i = 0; i < uptakeRates.length; i++ )
        {
            uptakeRates[i] = 0.0;
        }
    }

    public void scaleSecretion(double factor)
    {
        for( int i = 0; i < secretionRates.length; i++ )
        {
            secretionRates[i] *= factor;
            netExportRates[i] *= factor;
        }
    }

    public void scaleUptake(double factor)
    {
        for( int i = 0; i < uptakeRates.length; i++ )
        {
            uptakeRates[i] *= factor;
        }

    }

    public void setSecretionRate(int index, double val)
    {
        secretionRates[index] = val;
    }

    @Override
    public Secretion clone()
    {
        try
        {
            Secretion result = (Secretion)super.clone();
            result.secretionRates = this.secretionRates.clone();
            result.uptakeRates = this.uptakeRates.clone();
            result.saturationDensities = this.saturationDensities.clone();
            result.netExportRates = this.netExportRates.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }

    public String display()
    {
        List<String> secretion = new ArrayList<>();
        List<String> uptake = new ArrayList<>();
        List<String> export = new ArrayList<>();
        for( int i = 0; i < substrates.length; i++ )
        {
            if( secretionRates[i] != 0 )
                secretion.add( substrates[i] + ", rate " + secretionRates[i] );
            if( uptakeRates[i] != 0 )
                uptake.add( substrates[i] + ", rate " + uptakeRates[i] );
            if( netExportRates[i] != 0 )
                export.add( substrates[i] + ", rate " + netExportRates[i] );
        }
        if( secretion.isEmpty() && uptake.isEmpty() && export.isEmpty() )
            return "Secretion Disabled.\n--------------------------------";

        StringBuilder sb = new StringBuilder();
        sb.append( "Secretion:\n--------------------------------" );
        if( !secretion.isEmpty() )
            sb.append( "\n\tSecretes " + secretion.get( 0 ) );
        for( int i = 1; i < secretion.size(); i++ )
            sb.append( "\n\t         " + secretion.get( i ) );
        if( !uptake.isEmpty() )
            sb.append( "\n\tUptakes " + uptake.get( 0 ) );
        for( int i = 1; i < uptake.size(); i++ )
            sb.append( "\n\t        " + uptake.get( i ) );
        if( !export.isEmpty() )
            sb.append( "\n\tExports " + export.get( 0 ) );
        for( int i = 1; i < export.size(); i++ )
            sb.append( "\n\t        " + export.get( i ) );
        return sb.toString();
    }
}
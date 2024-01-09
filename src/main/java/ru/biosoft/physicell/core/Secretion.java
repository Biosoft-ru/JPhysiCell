package ru.biosoft.physicell.core;

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
    double[] secretionRates = new double[0];
    public double[] uptakeRates = new double[0];
    double[] saturationDensities = new double[0];
    double[] netExportRates = new double[0];

    public void sync(Microenvironment m)
    {
        int size = m.number_of_densities();
        secretionRates = VectorUtil.resize( secretionRates, size );//new double[m.number_of_densities()];//.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
        uptakeRates = VectorUtil.resize( uptakeRates, size );//.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
        saturationDensities = VectorUtil.resize( saturationDensities, size );//.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
        netExportRates = VectorUtil.resize( netExportRates, size );//.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
    }

    public void advance(BasicAgent cell, Phenotype phenotype, double dt)
    {
        // if this phenotype is not associated with a cell, exit 
        if( cell == null )
            return;

        // make sure the associated cell has the correct rate vectors 
        if( cell.secretionRates != secretionRates )
        {
            cell.secretionRates = secretionRates;
            cell.uptakeRates = uptakeRates;
            cell.saturationDensities = saturationDensities;
            cell.netExportRates = netExportRates;
            cell.setTotalVolume( phenotype.volume.total );
            cell.setUptakeConstants( dt );
        }

        // now, call the BioFVM secretion/uptake function 
        cell.simulateSecretionUptake( cell.getMicroenvironment(), dt );
    }

    public void set_all_secretion_to_zero()
    {
        for( int i = 0; i < secretionRates.length; i++ )
        {
            secretionRates[i] = 0.0;
            netExportRates[i] = 0.0;
        }
    }

    void set_all_uptake_to_zero()
    {
        for( int i = 0; i < uptakeRates.length; i++ )
        {
            uptakeRates[i] = 0.0;
        }
    }

    public void scale_all_secretion_by_factor(double factor)
    {
        for( int i = 0; i < secretionRates.length; i++ )
        {
            secretionRates[i] *= factor;
            netExportRates[i] *= factor;
        }
    }

    public void scale_all_uptake_by_factor(double factor)
    {
        for( int i = 0; i < uptakeRates.length; i++ )
        {
            uptakeRates[i] *= factor;
        }

    }

    // ease of access
    public void setSecretionRate(int index, double val)
    {
        //        int index = microenvironment.find_density_index( name );
        secretionRates[index] = val;
        //        return 0;
    }

    double uptake_rate(String name)
    {
        //            int index = microenvironment.find_density_index(name); 
        //            return uptake_rates[index]; 
        return 0;
    }

    double saturation_density(String name)
    {
        //            int index = microenvironment.find_density_index(name); 
        //            return saturation_densities[index];
        return 0;
    }

    double net_export_rate(String name)
    {
        //            int index = microenvironment.find_density_index(name); 
        //            return net_export_rates[index];
        return 0;
    }

    @Override
    public Secretion clone()
    {
        try
        {
            Secretion result = (Secretion)super.clone();
            result.secretionRates = this.secretionRates.clone();
            result.uptakeRates = this.secretionRates.clone();
            result.saturationDensities = this.saturationDensities.clone();
            result.netExportRates = this.netExportRates.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
}
package ru.biosoft.biofvm.cell;

import ru.biosoft.biofvm.BasicAgent;
import ru.biosoft.biofvm.Microenvironment;

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
    Microenvironment pMicroenvironment;

    double[] secretion_rates;
    public double[] uptake_rates;
    double[] saturation_densities;
    double[] net_export_rates;

    public Secretion()
    {
        pMicroenvironment = Microenvironment.get_default_microenvironment();

        sync_to_current_microenvironment();
    }

    public void sync_to_current_microenvironment()
    {
        if( pMicroenvironment != null )
        {
            sync_to_microenvironment( pMicroenvironment );
        }
        else
        {
            secretion_rates = new double[0];
            uptake_rates = new double[0];
            saturation_densities = new double[0];
            net_export_rates = new double[0];
        }
    }

    public void sync_to_microenvironment(Microenvironment pNew_Microenvironment)
    {
        pMicroenvironment = pNew_Microenvironment;

        secretion_rates = new double[pMicroenvironment.number_of_densities()];//.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
        uptake_rates = new double[pMicroenvironment.number_of_densities()];//.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
        saturation_densities = new double[pMicroenvironment.number_of_densities()];//.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
        net_export_rates = new double[pMicroenvironment.number_of_densities()];//.resize( pMicroenvironment->number_of_densities() , 0.0 ); 
    }

    public void advance(BasicAgent pCell, Phenotype phenotype, double dt)
    {
        // if this phenotype is not associated with a cell, exit 
        if( pCell == null )
            return;

        // if there is no microenvironment, attempt to sync. 
        if( pMicroenvironment == null )
        {
            // first, try the cell's microenvironment
            if( pCell.getMicroenvironment() != null )
            {
                sync_to_microenvironment( pCell.getMicroenvironment() );
            }
            // otherwise, try the default microenvironment
            else
            {
                sync_to_microenvironment( Microenvironment.get_default_microenvironment() );
            }

            // if we've still failed, return. 
            if( pMicroenvironment == null )
                return;
        }

        // make sure the associated cell has the correct rate vectors 
        if( pCell.secretionRates != secretion_rates )
        {
            pCell.secretionRates = secretion_rates;
            pCell.uptakeRates = uptake_rates;
            pCell.saturationDensities = saturation_densities;
            pCell.netExportRates = net_export_rates;
            pCell.set_total_volume( phenotype.volume.total );
            pCell.setUptakeConstants( dt );
        }

        // now, call the BioFVM secretion/uptake function 
        pCell.simulateSecretionUptake( pMicroenvironment, dt );
    }

    public void set_all_secretion_to_zero()
    {
        for( int i = 0; i < secretion_rates.length; i++ )
        {
            secretion_rates[i] = 0.0;
            net_export_rates[i] = 0.0;
        }
    }

    void set_all_uptake_to_zero()
    {
        for( int i = 0; i < uptake_rates.length; i++ )
        {
            uptake_rates[i] = 0.0;
        }
    }

    public void scale_all_secretion_by_factor(double factor)
    {
        for( int i = 0; i < secretion_rates.length; i++ )
        {
            secretion_rates[i] *= factor;
            net_export_rates[i] *= factor;
        }
    }

    public void scale_all_uptake_by_factor(double factor)
    {
        for( int i = 0; i < uptake_rates.length; i++ )
        {
            uptake_rates[i] *= factor;
        }

    }

    // ease of access
    public void setSecretionRate(int index, double val)
    {
        //        int index = microenvironment.find_density_index( name );
        secretion_rates[index] = val;
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
            result.secretion_rates = this.secretion_rates.clone();
            result.uptake_rates = this.secretion_rates.clone();
            result.saturation_densities = this.saturation_densities.clone();
            result.net_export_rates = this.net_export_rates.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
}
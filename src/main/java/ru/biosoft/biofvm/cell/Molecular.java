package ru.biosoft.biofvm.cell;

import ru.biosoft.biofvm.BasicAgent;
import ru.biosoft.biofvm.Microenvironment;
import ru.biosoft.biofvm.VectorUtil;

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
public class Molecular implements Cloneable
{
    Microenvironment pMicroenvironment;


    // we'll set this to replace BioFVM's version       
    double[] internalized_total_substrates;

    // for each substrate, a fraction 0 <= f <= 1 of the 
    // total internalized substrate is released back inot
    // the environment at death 
    double[] fraction_released_at_death;

    // for each substrate, a fraction 0 <= f <= 1 of the 
    // total internalized substrate is transferred to the  
    // predatory cell when ingested 
    double[] fraction_transferred_when_ingested;

    public Molecular()
    {
        pMicroenvironment = Microenvironment.get_default_microenvironment();
        sync_to_current_microenvironment();
    }

    void sync_to_current_microenvironment()
    {
        if( pMicroenvironment != null )
        {
            sync_to_microenvironment( pMicroenvironment );
        }
        else
        {
            internalized_total_substrates = new double[0];//.resize( 0, 0.0 );
            fraction_released_at_death = new double[0];//.resize( 0, 0.0 );
            fraction_transferred_when_ingested = new double[0];//.resize( 0, 0.0 );
        }
        return;
    }

    void sync_to_microenvironment(Microenvironment pNew_Microenvironment)
    {
        pMicroenvironment = pNew_Microenvironment;
        int number_of_densities = pMicroenvironment.number_of_densities();
        //        internalized_total_substrates.resize( number_of_densities, 0.0 );
        //        fraction_released_at_death.resize( number_of_densities, 0.0 );
        //        fraction_transferred_when_ingested.resize( number_of_densities, 0.0 );
        internalized_total_substrates = VectorUtil.resize( internalized_total_substrates, number_of_densities );
        fraction_released_at_death = VectorUtil.resize( internalized_total_substrates, number_of_densities );
        fraction_transferred_when_ingested = VectorUtil.resize( fraction_transferred_when_ingested, number_of_densities );
    }

    void sync_to_cell(BasicAgent pCell)
    {
        //                delete pCell.internalized_substrates;
        //                pCell.internalized_substrates = &internalized_total_substrates;
        pCell.internalizedSubstrates = internalized_total_substrates;

        //                delete pCell.fraction_released_at_death;
        //                pCell.fraction_released_at_death = &fraction_released_at_death; 
        pCell.fraction_released_at_death = fraction_released_at_death;

        //                delete pCell.fraction_transferred_when_ingested; 
        //                pCell.fraction_transferred_when_ingested = &fraction_transferred_when_ingested; 
        pCell.fraction_transferred_when_ingested = fraction_transferred_when_ingested;
    }

    // ease of access 
    //            double& internalized_total_substrate( String name )
    //            {
    //                int index = microenvironment.find_density_index(name); 
    //                return internalized_total_substrates[index]; 
    //            }

    /*
    void Molecular::advance( Basic_Agent* pCell, Phenotype& phenotype , double dt )
    {
        // if this phenotype is not associated with a cell, exit 
        if( pCell == NULL )
        { return; }
    
        // if there is no microenvironment, attempt to sync. 
        if( pMicroenvironment == NULL )
        {
            // first, try the cell's microenvironment
            if( pCell->get_microenvironment() )
            {
                sync_to_microenvironment( pCell->get_microenvironment() ); 
            }
            // otherwise, try the default microenvironment
            else
            {
                sync_to_microenvironment( get_default_microenvironment() ); 
            }
    
            // if we've still failed, return. 
            if( pMicroenvironment == NULL ) 
            {
                return; 
            }
        }
    
        // make sure the associated cell has the correct rate vectors 
        if( pCell->internalized_substrates != &internalized_substrates )
        {
            // copy the data over 
            internalized_substrates = *(pCell->internalized_substrates);
            // remove the BioFVM copy 
            delete pCell->internalized_substrates; 
            // point BioFVM to this one  
            pCell->internalized_substrates = &internalized_substrates; 
        }
    
        // now, call the functions 
    //  if( pCell->functions.internal_substrate_function )
    //  { pCell->functions.internal_substrate_function( pCell,phenotype,dt);  }
    //  if( pCell->functions.molecular_model_function )
    //  { pCell->functions.molecular_model_function( pCell,phenotype,dt);  }
    
    
        return; 
    }
    */
    @Override
    public Molecular clone()
    {
        try
        {
            Molecular result = (Molecular)super.clone();
            result.internalized_total_substrates = this.internalized_total_substrates.clone();
            result.fraction_released_at_death = this.fraction_released_at_death.clone();
            result.fraction_transferred_when_ingested = this.fraction_transferred_when_ingested.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }
}
package ru.biosoft.physicell.core;

import java.util.ArrayList;
import java.util.List;

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
public class CellInteractions implements Cloneable
{
    public double apoptotic_phagocytosis_rate;
    public double necrotic_phagocytosis_rate;
    public double other_dead_phagocytosis_rate;

    public Cell pAttackTarget = null; 
    public double total_damage_delivered = 0.0; 
    
    @Deprecated
    public double deadPhagocytosisRate; // phagocytosis parameters (e.g., macrophages)
    public double[] livePhagocytosisRates; // attack parameters (e.g., T cells)
    public double[] attackRates; // do I attack cell type j? 
    public double[] immunogenicities; // how immnogenic am I to cell type j?  
    public double damageRate;
    public double attack_duration = 30; // a typical attack duration for a T cell using perforin/granzyme is ~30 minutes
    public double[] fusionRates; // cell fusion parameters 
    private String[] names;

    public CellInteractions()
    {
        deadPhagocytosisRate = 0.0;
        livePhagocytosisRates = new double[] {0.0};
        damageRate = 1.0;
        attackRates = new double[] {0.0};
        immunogenicities = new double[] {1};
        fusionRates = new double[] {0.0};
    }

    public void initialize(Model model)
    {

        int size = model.getDefinitionsCount();
        names = new String[size];
        for( int i = 0; i < size; i++ )
            names[i] = model.getCellDefinition( i ).name;
        livePhagocytosisRates = VectorUtil.resize( livePhagocytosisRates, size );
        attackRates = VectorUtil.resize( attackRates, size );
        fusionRates = VectorUtil.resize( fusionRates, size );
        immunogenicities = VectorUtil.resize( immunogenicities, size, 1 );
    }

    public double getLivePhagocytosisRate(int type)
    {
        return livePhagocytosisRates[type];
    }

    public double getAttackRate(int type)
    {
        return attackRates[type];
    }

    public double getImmunogenicity(int type)
    {
        return immunogenicities[type];
    }

    public double getFusionRate(int type)
    {
        return fusionRates[type];
    }

    @Override
    public CellInteractions clone()
    {
        try
        {
            CellInteractions result = (CellInteractions)super.clone();
            result.livePhagocytosisRates = this.livePhagocytosisRates.clone();
            result.attackRates = this.attackRates.clone();
            result.immunogenicities = this.immunogenicities.clone();
            result.fusionRates = this.fusionRates.clone();
            return result;
        }
        catch( CloneNotSupportedException e )
        {
            throw ( new InternalError( e ) );
        }
    }

    public String display()
    {
        StringBuilder sb = new StringBuilder();
        List<String> attacks = new ArrayList<>();
        List<String> fusions = new ArrayList<>();
        List<String> phogocyte = new ArrayList<>();

        for( int i = 0; i < attackRates.length; i++ )
        {
            if( attackRates[i] != 0 )
                attacks.add( names[i] );
            if( fusionRates[i] != 0 )
                fusions.add( names[i] );
            if( livePhagocytosisRates[i] != 0 )
                phogocyte.add( names[i] );
        }

        if( attacks.isEmpty() && fusions.isEmpty() && phogocyte.isEmpty() )
        {
            sb.append( "Interactions Disabled." );
            sb.append( "\n--------------------------------" );
            return sb.toString();

        }

        sb.append( "Interactions:" );
        sb.append( "\n--------------------------------" );
        if( !attacks.isEmpty() )
            sb.append( "\n\tAttacks " );

        for( int i = 0; i < attacks.size(); i++ )
        {
            if( i > 0 )
                sb.append( "\n\t        " );
            sb.append( attacks.get( i ) );
        }

        if( !fusions.isEmpty() )
            sb.append( "\n\tFuses " );
        for( int i = 0; i < fusions.size(); i++ )
        {
            if( i > 0 )
                sb.append( "\n\t      " );
            sb.append( fusions.get( i ) );
        }

        if( !phogocyte.isEmpty() )
            sb.append( "\n\tIngests " );
        for( int i = 0; i < phogocyte.size(); i++ )
        {
            if( i > 0 )
                sb.append( "\n\t         " );
            sb.append( phogocyte.get( i ) );
        }
        return sb.toString();
    }
}
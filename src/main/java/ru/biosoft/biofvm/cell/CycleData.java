package ru.biosoft.biofvm.cell;

import java.util.ArrayList;
import java.util.List;

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
public class CycleData
{
    private CycleModel cycleModel;

    String timeUnits;

    public List<List<Double>> transitionRates;

    public List<List<Double>> basicRates;

    public int currentPhaseIndex;
    public double elapsedTimePhase;

    public CycleData(CycleModel model)
    {
        cycleModel = model;
        timeUnits = "min";
        transitionRates = new ArrayList<>();
        basicRates = new ArrayList<>();
        currentPhaseIndex = 0;
        elapsedTimePhase = 0.0;
    }

    private void resizeListList(List<List<Double>> list, int n)
    {
        if( n < list.size() )
        {
            for( int i = n; i < list.size(); i++ )
                list.remove( i );
        }
        else if( n > list.size() )
        {
            for( int i = 0; i < n - list.size(); i++ )
                list.add( new ArrayList<Double>() );
        }
    }

    private void resizeList(List<Double> list, int n)
    {
        if( n < list.size() )
        {
            for( int i = n; i < list.size(); i++ )
                list.remove( i );
        }
        else if( n > list.size() )
        {
            for( int i = 0; i < n - list.size(); i++ )
                list.add( 0.0 );
        }
    }

    public void sync_to_cycle_model()
    {
        int n = cycleModel.phases.size();
        resizeListList( transitionRates, n );
        resizeListList( basicRates, n );
        for( int i = 0; i < cycleModel.phase_links.size(); i++ )
        {
            int size = cycleModel.phase_links.get( i ).size();
            for( int j = 0; j < size; j++ )
            {
                List<Double> rates = transitionRates.get( i );
                resizeList( rates, size );
                resizeList( basicRates.get( i ), size );
            }
        }
    }

    double getTransitionRate(int start_phase_index, int end_phase_index)
    {
        int index = cycleModel.inverse_index_maps.get( start_phase_index ).get( end_phase_index );
        return transitionRates.get( start_phase_index ).get( index );
    }

    public void setTransitionRate(int start_phase_index, int end_phase_index, double rate)
    {
        int index = cycleModel.inverse_index_maps.get( start_phase_index ).get( end_phase_index );
        transitionRates.get( start_phase_index ).set( index, rate );
        basicRates.get( start_phase_index ).set( index, rate );
    }

    public void modifyTransitionRate(int start_phase_index, int end_phase_index, double multiplier)
    {
        int index = cycleModel.inverse_index_maps.get( start_phase_index ).get( end_phase_index );
        double basic = this.basicRates.get( start_phase_index ).get( index );
        transitionRates.get( start_phase_index ).set( index, basic * multiplier );
    }

    double getExitRate(int phase_index)
    {
        return transitionRates.get( phase_index ).get( 0 );
    }
    public void setExitRate(int phase_exit, double rate)
    {
        transitionRates.get( phase_exit ).set( 0, rate );
    }

    public Phase currentPhase()
    {
        return cycleModel.phases.get( currentPhaseIndex );
    }

    public CycleData clone(CycleModel model)
    {
        CycleData result = new CycleData( model );
        result.timeUnits = this.timeUnits;
        result.currentPhaseIndex = this.currentPhaseIndex;
        result.elapsedTimePhase = this.elapsedTimePhase;
        result.transitionRates = new ArrayList<>();
        result.basicRates = new ArrayList<>();
        for( int i = 0; i < transitionRates.size(); i++ )
        {
            result.transitionRates.add( new ArrayList<Double>( transitionRates.get( i ) ) );
            result.basicRates.add( new ArrayList<Double>( transitionRates.get( i ) ) );
        }
        return result;
    }
}
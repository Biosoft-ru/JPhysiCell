package ru.biosoft.physicell.fba;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;

import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;

import biouml.plugins.fbc.ApacheModel;

public class FBA_model
{
    public ApacheModel fbcModel;

    private Map<String, Integer> reactionIndex;
    private Map<String, LinearConstraint> lowerBoundIndex;

    public FBA_model()
    {
    }

    public double getFlux(String reaction)
    {
        return fbcModel.getOptimValue( reaction );
    }

    public void setApacheModel(ApacheModel model)
    {
        this.fbcModel = model;
        initReactionMap( model );
        initLowerBoundMap( model );
    }

    public void initReactionMap(ApacheModel model)
    {
        reactionIndex = new HashMap<>();
        String[] names = model.getReactionNames();
        for( int i = 0; i < names.length; i++ )
            reactionIndex.put( names[i], i );
    }

    public void initLowerBoundMap(ApacheModel model)
    {
        lowerBoundIndex = new HashMap<>();
        Collection<LinearConstraint> constraints = model.getConstraints();
        for( String reaction : fbcModel.getReactionNames() )
        {
            int index = reactionIndex.get( reaction );
            for( LinearConstraint lc : constraints )
            {
                if( isLowerBound( index, lc ) )
                {
                    lowerBoundIndex.put( reaction, lc );
                    break;
                }
            }
        }
    }

    protected boolean isLowerBound(int reactionIndex, LinearConstraint constraint)
    {
        Relationship relation = constraint.getRelationship();
        if( !relation.equals( Relationship.GEQ ) )
            return false;
        RealVector rv = constraint.getCoefficients();
        if( rv.getEntry( reactionIndex ) != 1 )
            return false;
        for( int i = 0; i < rv.getMaxIndex(); i++ )
        {
            if( i != reactionIndex && rv.getEntry( i ) != 0 )
                return false;
        }
        return true;
    }

    public void setReactionLowerBound(String rId, double lowerBound)
    {
        LinearConstraint lc = lowerBoundIndex.get( rId );
        LinearConstraint newLc = new LinearConstraint( lc.getCoefficients(), Relationship.GEQ, lowerBound );
        fbcModel.getConstraints().remove( lc );
        fbcModel.getConstraints().add( newLc );
    }

    public boolean getSolutionStatus()
    {
        return true; //@TODO
    }

    public double getObjectiveValue()
    {
        return fbcModel.getValueObjFunc();
    }

    public void runFBA()
    {
        fbcModel.setLogLevel( Level.SEVERE );
        fbcModel.optimize();
        //        System.out.println( "Solved " + fbcModel.getValueObjFunc() );
        //@TODO
        //        std::cout << "Running FBA... ";
        //        this.lp_model.primal();
        //        if( lp_model.isProvenOptimal() )
        //        {
        //            double[] columnPrimal = this.lp_model.primalColumnSolution();
        //            //            std::cout << "Optimal solution found!" << std::endl;
        //            for( FBA_reaction reaction : this.reactions )
        //            {
        //                int idx = this.reactionsIndexer.get( reaction.getId() );
        //                double v = columnPrimal[idx];
        //                reaction.setFluxValue( v );
        //            }
        //        }
        //        else
        //        {
        //            for( FBA_reaction reaction : this.reactions )
        //            {
        //                reaction.setFluxValue( 0.0 );
        //            }
        //            //            std::cout << "Primal infeasible" << std::endl;
        //        }

    }

    public FBA_model clone()
    {
        FBA_model result = new FBA_model();
        result.setApacheModel( (ApacheModel)fbcModel.clone() );
        return result;
    }

    //
    //    void initLpModel()
    //    {
    //
    //        //        int n_rows = this.getNumMetabolites();
    //        //        int n_cols = this.getNumReactions();
    //        //
    //        //        //this.lp_model = new ClpSimplex();
    //        //        handler = new CoinMessageHandler(null);
    //        //        
    //        //        handler.setLogLevel(0);
    //        //        lp_model.passInMessageHandler(handler);
    //        //
    //        //        CoinPackedMatrix matrix;
    //        //        matrix.setDimensions(n_rows, 0);
    //        //
    //        //        double[] row_lb = new double[n_rows]; //the row lower bounds
    //        //        double[] row_ub = new double[n_rows]; //the row upper bounds
    //        //        double[] col_lb = new double[n_cols]; //the column lower bounds
    //        //        double[] col_ub = new double[n_cols]; //the column upper bounds
    //        //        double[] objective = new double[n_cols]; //the objective coefficients
    //        //
    //        //        for(int i=0; i< n_rows; i++)
    //        //        {
    //        //            row_lb[i] = 0;
    //        //            row_ub[i] = 0;
    //        //        }
    //        //        for(FBA_reaction reaction: this.reactions)
    //        //        {
    //        //            int col_idx = this.reactionsIndexer.get( reaction.getId() );
    //        //            col_lb[col_idx] = reaction.getLowerBound();
    //        //            col_ub[col_idx] = reaction.getUpperBound();
    //        //            objective[col_idx] = reaction.getObjectiveCoefficient();
    //        //
    //        //            Map<FBA_metabolite, Double> metabolites = reaction.getMetabolites();
    //        //            CoinPackedVector col;
    //        //            //            for(auto it=metabolites.begin(); it!=metabolites.end(); it++)
    //        //            //            {
    //        //            //                const FBA_metabolite* metabolite = it.first;
    //        //            //                double stoich_coeff = it.second;
    //        //            for( Entry<FBA_metabolite, Double> entry : metabolites.entrySet() )
    //        //            {
    //        //                FBA_metabolite metabolite = entry.getKey();
    //        //                Double stoich_coeff = entry.getValue();
    //        //                int row_idx = this.metaboliteIndexer.get( metabolite.getId() );
    //        //                col.insert(row_idx, stoich_coeff);
    //        //            }
    //        //            matrix.appendCol(col);
    //        //        }
    //        //
    //        //        this.lp_model.loadProblem(matrix, col_lb, col_ub, objective, row_lb, row_ub);
    //        //        this.lp_model.setOptimizationDirection(-1);
    //        //
    //        //        //        delete col_lb;
    //        //        //        delete col_ub;
    //        //        //        delete row_lb;
    //        //        //        delete row_ub;
    //        //        //        delete objective;
    //        //
    //        //        this.is_initialized = true;
    //    }
}
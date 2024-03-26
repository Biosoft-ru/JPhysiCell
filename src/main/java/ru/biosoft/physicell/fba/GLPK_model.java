package ru.biosoft.physicell.fba;

import java.util.logging.Level;

import org.gnu.glpk.GLPK;
import org.gnu.glpk.GLPKConstants;

import biouml.plugins.fbc.GLPKModel;

public class GLPK_model
{
    public GLPKModel fbcModel;

    public GLPK_model()
    {
    }

    public double getFlux(String reaction)
    {
        return fbcModel.getOptimValue( reaction );
    }

    public void setGLPKModel(GLPKModel model)
    {
        this.fbcModel = model;
        model.setLogLevel( Level.SEVERE );
        GLPK.glp_term_out( GLPKConstants.GLP_OFF );
    }

    public void setReactionLowerBound(String rId, double lowerBound)
    {
        fbcModel.setLowerBound( lowerBound, rId );
    }

    public boolean getSolutionStatus()
    {
        int status = GLPK.glp_get_prim_stat( fbcModel.problem );
        return status == GLPK.GLP_FEAS;
    }

    public double getObjectiveValue()
    {
        return fbcModel.getValueObjFunc();
    }

    public void runFBA()
    {
        fbcModel.optimize();
    }

    public GLPK_model clone()
    {
        GLPK_model result = new GLPK_model();
        result.setGLPKModel( (GLPKModel)fbcModel.clone() );
        return result;
    }
}
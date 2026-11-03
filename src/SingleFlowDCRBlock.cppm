/*--------------------------------------------------------------------------*/
/*---------------------- File SingleFlowDCRBlock.cppm ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.singleflowdcrblock, which exports every
 * name that the headers of SingleFlowDCRBlock declare at namespace scope,
 * together with the modules it depends on; a translation unit can then write
 *
 *     import smspp.singleflowdcrblock;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Laura Galli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Luca Mencarelli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Sorbera \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Laura Galli, Luca Mencarelli,
 *                      Enrico Sorbera
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "DCR_SPT.h"
#include "MultiFlowDCRBlock.h"
#include "SingleFlowDCRBendersSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.singleflowdcrblock;

export import smspp;

export namespace SMSpp_di_unipi_it {

 // MultiFlowDCRBlock.h
 using SMSpp_di_unipi_it::MultiFlowDCRBlock;
 using SMSpp_di_unipi_it::MultiSubset;
 using SMSpp_di_unipi_it::MultiVector;
 using SMSpp_di_unipi_it::Vec_Bool;
 using SMSpp_di_unipi_it::Vec_double;
 using SMSpp_di_unipi_it::c_RHSValue;
 using SMSpp_di_unipi_it::c_Vec_double;

 // SingleFlowDCRBendersSolver.h
 using SMSpp_di_unipi_it::SingleFlowDCRBendersSolver;

 // SingleFlowDCRBlock.h
 using SMSpp_di_unipi_it::DCRSolution;
 using SMSpp_di_unipi_it::SingleFlowDCRBlock;
 using SMSpp_di_unipi_it::SingleFlowDCRBlockMod;
 using SMSpp_di_unipi_it::SingleFlowDCRBlockRngdMod;
 using SMSpp_di_unipi_it::SingleFlowDCRBlockSbstMod;
 using SMSpp_di_unipi_it::Vec_SingleFlowDCRBlock;
 using SMSpp_di_unipi_it::Vec_SingleFlowDCRBlock_it;
 using SMSpp_di_unipi_it::c_Vec_SingleFlowDCRBlock;
 using SMSpp_di_unipi_it::c_Vec_SingleFlowDCRBlock_it;
 using SMSpp_di_unipi_it::p_SingleFlowDCRBlock;
}

/*--------------------------------------------------------------------------*/
/*-------------------- End File SingleFlowDCRBlock.cppm --------------------*/
/*--------------------------------------------------------------------------*/

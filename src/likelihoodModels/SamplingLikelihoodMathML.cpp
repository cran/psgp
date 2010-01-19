#include "SamplingLikelihoodMathML.h"

#include <cmath>

using namespace std;
using namespace itpp;
XERCES_CPP_NAMESPACE_USE

SamplingLikelihoodMathML::SamplingLikelihoodMathML(string MathMLString)
{
	bool errorOccurred = false;

	try
	{
		XMLPlatformUtils::Initialize();
	}
	catch (const XMLException& toCatch)
	{
		cerr << "Error during initialization! :\n" << XMLString::transcode(toCatch.getMessage()) << endl;
		return;
	}

	static const XMLCh gLS[] = { chLatin_L, chLatin_S, chNull };
	DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(gLS);
	DOMLSParser       *parser = ((DOMImplementationLS*)impl)->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
	DOMConfiguration  *config = parser->getDomConfig();

	config->setParameter(XMLUni::fgDOMNamespaces, false);
	config->setParameter(XMLUni::fgXercesSchema, false);
	config->setParameter(XMLUni::fgXercesSchemaFullChecking, false);
	config->setParameter(XMLUni::fgDOMValidateIfSchema, true);
	config->setParameter(XMLUni::fgDOMDatatypeNormalization, true);

	const char* xmlMessageContent = MathMLString.c_str();
	MemBufInputSource* xmlMessage = new MemBufInputSource((const XMLByte*)xmlMessageContent, strlen(xmlMessageContent), "MathMLString", false);
	Wrapper4InputSource *wrapper=new Wrapper4InputSource(xmlMessage,false);

	XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument *doc;

	try
	{
		parser->resetDocumentPool();
		doc = parser->parse((DOMLSInput *)wrapper);

		if(doc)
	   {
			MathMLParser mp;
			expressionTree = mp.parseExpression((DOMNode *)doc->getDocumentElement());
			cout << "Parsed Expression: ";
			cout << expressionTree->toString() << endl;
		}
	}
	catch (const XMLException& toCatch)
	{
		cerr << "\nError during parsing: Exception message is:  \n" << XMLString::transcode(toCatch.getMessage()) << "\n" << endl;
		errorOccurred = true;
	}
	catch (const DOMException& toCatch)
	{
		const unsigned int maxChars = 2047;
		XMLCh errText[maxChars + 1];
		cerr << "\nDOM Error during parsing: DOMException code is:  " << toCatch.code << endl;
		if (DOMImplementation::loadDOMExceptionMsg(toCatch.code, errText, maxChars))
		{
			cerr << "Message is: " << XMLString::transcode(errText) << endl;
		}
		errorOccurred = true;
	}
	catch (...)
	{
		cerr << "\nUnexpected exception during parsing\n";
		errorOccurred = true;
	}

	parser->release();
	XMLPlatformUtils::Terminate();
}


SamplingLikelihoodMathML::SamplingLikelihoodMathML(DOMDocument *doc)
{
	MathMLParser mp;
	expressionTree = mp.parseExpression((DOMNode *)doc->getDocumentElement());
	cout << "Parsed Expression: ";
	cout << expressionTree->toString() << endl;
}

SamplingLikelihoodMathML::SamplingLikelihoodMathML(ExpressionNode *tree)
{
	expressionTree = tree;
	cout << "Parsed Expression: ";
	cout << expressionTree->toString() << endl;
}

SamplingLikelihoodMathML::~SamplingLikelihoodMathML()
{
	delete expressionTree;
}

double SamplingLikelihoodMathML::modelFunction(const double x) const
{
	return expressionTree->evaluate(x);
}

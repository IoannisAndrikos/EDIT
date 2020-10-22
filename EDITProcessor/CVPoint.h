//Copyright (c)2012 FORTH  Greece
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE

#pragma once
#include<math.h>
using namespace System::Collections::Generic;

using namespace System::Xml::Serialization;
namespace EDITCore
{
	public  ref class CVPoint
	{

	private:
		double x, y;// z;

	public:

		double GetX() { return x; }
		double GetY() { return y; }
		//double GetZ() { return z; }

		void SetX(double x1) { x = x1; }
		void SetY(double y1) { y = y1; }
		//void SetZ(double z1) { z = z1;; }

		CVPoint(double x1, double y1) { x = x1; y = y1; } // , double z1 //  z = z1;


		CVPoint(void)
		{

		}

		double Dist(CVPoint ^p, double sx, double sy) //, double sz
		{


			return sqrt((p->GetX() - x)*(p->GetX() - x)*sx*sx + (p->GetY() - y)*(p->GetY() - y)*sy*sy);// +(p->GetZ() - z)*(p->GetZ() - z)*sz*sz);

		}


		[XmlAttributeAttribute]
		property System::Double X
		{
			double get()
			{
				return x;
			}
			void set(double x1)
			{
				x = x1;
			}
		}

		[XmlAttributeAttribute]
		property System::Double Y
		{
			double get()
			{
				return y;
			}
			void set(double y1)
			{
				y = y1;
			}

		}

		/*[XmlAttributeAttribute]
		property System::Double Z
		{
			double get()
			{
				return z;
			}
			void set(double z1)
			{
				z = z1;
			}

		}*/
	};

	public ref class EDITResponse
	{
	private:
		static EDITResponse instance;
		EDITResponse() {}
		EDITResponse(const EDITResponse%) { throw gcnew System::InvalidOperationException("singleton cannot be copy-constructed"); }

		System::String^ errorMessage = "";
		List<List<EDITCore::CVPoint^>^>^ allFramesData = gcnew List<List<EDITCore::CVPoint^>^>();
		List<EDITCore::CVPoint^>^ uniqueFrameData = gcnew List<EDITCore::CVPoint^>();
		List<double>^ numericData = gcnew List<double>();
		System::String^ path = "";
		List <System::String^>^ paths = gcnew List <System::String^>();

	public:
		static property EDITResponse^ Instance { EDITResponse^ get() { return % instance; } }

		bool isSuccessful() {
			if (this->errorMessage == "success") {
				return true;
			}
			else {
				return false;
			}
		}

		void setSuccessOrFailure(System::String^ errorMessage) {
			this->allFramesData->Clear();
			this->uniqueFrameData->Clear();
			this->numericData->Clear();
			this->paths->Clear();
			this->errorMessage = errorMessage;
		}

		System::String^ getFailure() {
			return this->errorMessage;
		}

		void setData(List<List<EDITCore::CVPoint^>^>^ data) {
			this->allFramesData = data;
		}

		void setData(List<EDITCore::CVPoint^>^ data) {
			this->uniqueFrameData = data;
		}

		void setData(List<double>^ data) {
			this->numericData = data;
		}

		void setData(System::String^ path) {
			this->path = path;
		}

		void setData(List <System::String^>^ paths) {
			this->paths = paths;
		}

		List<List<EDITCore::CVPoint^>^>^ getAllFramesData() {
			return this->allFramesData;
		}

		List<EDITCore::CVPoint^>^ getUniqueFramesData() {
			return this->uniqueFrameData;
		}

		List<double>^ getNumericData() {
			return this->numericData;
		}

		System::String^ getPath() {
			return this->path;
		}

		List <System::String^>^ getPaths() {
			return this->paths;
		}
			
	};
}
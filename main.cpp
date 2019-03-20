#include "pch.h"
#include <iostream>
#include "mfe_nonlinear.h"

int main()
{
	NonlinearTask T;
	T.init();
	T.setParams();
	T.solve();

    std::cout << "Do it again!";
}


#include "testMainFile_StateDerivativeFunctionRKF.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    testMainFile_StateDerivativeFunctionRKF w;
    w.show();

    return a.exec();
}

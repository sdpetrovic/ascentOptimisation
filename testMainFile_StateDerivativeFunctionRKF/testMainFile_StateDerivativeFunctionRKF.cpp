#include "testMainFile_StateDerivativeFunctionRKF.h"
#include "ui_testMainFile_StateDerivativeFunctionRKF.h"

testMainFile_StateDerivativeFunctionRKF::testMainFile_StateDerivativeFunctionRKF(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::testMainFile_StateDerivativeFunctionRKF)
{
    ui->setupUi(this);
}

testMainFile_StateDerivativeFunctionRKF::~testMainFile_StateDerivativeFunctionRKF()
{
    delete ui;
}

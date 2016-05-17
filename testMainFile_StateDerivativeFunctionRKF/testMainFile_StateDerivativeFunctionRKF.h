#ifndef TESTMAINFILE_STATEDERIVATIVEFUNCTIONRKF_H
#define TESTMAINFILE_STATEDERIVATIVEFUNCTIONRKF_H

#include <QMainWindow>

namespace Ui {
class testMainFile_StateDerivativeFunctionRKF;
}

class testMainFile_StateDerivativeFunctionRKF : public QMainWindow
{
    Q_OBJECT

public:
    explicit testMainFile_StateDerivativeFunctionRKF(QWidget *parent = 0);
    ~testMainFile_StateDerivativeFunctionRKF();

private:
    Ui::testMainFile_StateDerivativeFunctionRKF *ui;
};

#endif // TESTMAINFILE_STATEDERIVATIVEFUNCTIONRKF_H

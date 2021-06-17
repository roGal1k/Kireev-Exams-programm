#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include <QVector>
#include <QLabel>
#include <QLineEdit>
#include <QGraphicsView>
#include <QGraphicsTextItem>

QT_BEGIN_NAMESPACE
namespace Ui { class Widget; }
QT_END_NAMESPACE

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();

private slots:
    void on_pushButton_clicked();
    void table(double &W_graph, double &T_graph);

private:
    Ui::Widget *ui;
    QGraphicsScene *scene;
    QGraphicsTextItem *text;
    QVector <double> *x;
    QVector <double> *y;
    //double x=0;
    //double y=0;
};
#endif // WIDGET_H

from django.urls import path
from .views import RetrieveDBSum, RetrieveDBData, RetrieveExtData , RetrieveDataSumByUser, DataUpload, Download, Update, Remove

urlpatterns = [
    path('retrieve/', RetrieveDBData.as_view()),
    path('retrieve-ext/', RetrieveExtData.as_view()),
    path('retrieve-by-user/', RetrieveDataSumByUser.as_view()),
    path('retrieve-db-sum/', RetrieveDBSum.as_view()),
    path('upload/', DataUpload.as_view()),
    path('download/', Download.as_view()),
    path('update/', Update.as_view()),
    path('remove/', Remove.as_view())
]



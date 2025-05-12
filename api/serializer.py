from rest_framework import serializers
from rest_framework_simplejwt.serializers import TokenObtainPairSerializer
from .models import CustomUser
from django.contrib.auth import authenticate
from django.utils.translation import gettext_lazy as _
from .models import Variant, Frequency, Submission,CNV, VarCounts

class VarSerializer(serializers.ModelSerializer):
    variation_id = serializers.CharField(required = False)
    chromosome = serializers.CharField(required = False)
    position = serializers.IntegerField(required = False)
    ref_allele = serializers.CharField(required = False)
    alt_allele = serializers.CharField(required = False)
    gene_name = serializers.CharField(required = False)
    consequence = serializers.CharField(required=False, allow_null = True)
    submission_id = serializers.IntegerField(required=False)

    class Meta:
        model = Variant
        fields = ('variation_id','chromosome','position','ref_allele','alt_allele','gene_name','consequence','submission_id')

class CNVarSerializer(serializers.ModelSerializer):

    submission_id = serializers.IntegerField()

    class Meta:
        model = CNV
        fields = ('variation_id','chromosome','start_pos', 'end_pos', 'site_count',
                  'consequence','submission_id','last_updated')

class FreqSerializer(serializers.ModelSerializer):

    variation_id = serializers.CharField()

    class Meta:
        model = Frequency
        fields = ('variation_id','homo_count','het_count','last_updated')

class SubSerializer(serializers.ModelSerializer):

    class Meta:
        model = Submission
        fields = ('submission_id','no_individuals','submission_date','username_id')

class UserSerializer(serializers.ModelSerializer):
    
    class Meta:
        model = CustomUser

        fields = ['password','last_login','is_superuser',
                  'username','first_name','last_name','is_staff',
                  'is_active','date_joined','email']

class CustomTokenObtainPairSerializer(TokenObtainPairSerializer):

    def validate(self, attrs):

        email = attrs.get('email')
        password = attrs.get('password')

        if email and password:
            user = authenticate(request=self.context.get("request"), email=email, password=password)

            if not user:
                raise serializers.ValidationError(_("No active account with given credentials"), code = "authorization")
        
        else:
            raise serializers.ValidationError(_("Must include 'email' and 'password'."), code = 'authorization')
        
        serialized_user = UserSerializer(user).data 
        token = self.get_token(user)

        return({
            "user":serialized_user,
            "refresh":str(token),
            "access":str(token.access_token)
        })

class VarCountSerializer(serializers.ModelSerializer):

    class Meta:

        model = VarCounts

        fields = ['entry_id', 'cnv_count','ssv_count','exo_count',
                  'intro_count','inter_count','reg_count','spli_count','intra_count']
